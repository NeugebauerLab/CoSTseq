import numpy as np
import pickle
from tqdm import tqdm
from src.bam_utils import reverse_complement
from src.dms_utils import gini, normalize_cov
from src.pro_utils import moving_average, c_freq, binned_signal
from scipy.stats import iqr, pearsonr
from Bio import SeqIO
import matplotlib.pyplot as plt
import matplotlib.colors as mplc
from matplotlib import cm
import pandas as pd
from copy import deepcopy
import gffutils
from collections import defaultdict
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

class DMS_Profile():
    """
    Class for handling DMS (Dimethyl Sulfate) profiling data.

    Parameters:
        - pkl_loc (str): Path to the pickle file containing DMS profile data.
        - sample_id (str): Identifier for the DMS profile sample.
        - genome (Genome): Genome object representing the genomic sequence.
        - min_cov (int): Minimum coverage threshold for filtering data (default: 700).
        - min_mut (int): Minimum mutation threshold for filtering data (default: 1).
        - norm_without_gu (bool): Option to normalize without considering 'G' and 'U' bases.
        - cov_corr_without_gu (bool): Option for coverage correction without considering 'G' and 'U' bases.
        - with_stats (bool): Option to load statistical information during initialization, requires a lot of memory.

    Methods:
        - __init__: Initialize a DMS_Profile object with specified parameters.
        - load_profile: Load DMS profile data from a pickle file.
        - get_mask: Create masks based on coverage and mutation thresholds.
        - calculate_rates: Calculate mutation rates with optional coverage correction.
        - plot_profile: Plot DMS profile for a specified genomic region.

    Attributes:
        - Various attributes for storing profile data and parameters.
    """

    def __init__(self, pkl_loc, sample_id, genome, min_cov=700, min_mut=1, norm_without_gu=False, cov_corr_without_gu=True, with_stats=False):
        # parameters
        self.min_cov = min_cov
        self.min_mut = min_mut

        self.pkl_loc = pkl_loc
        self.std_p   = None
        self.std_m   = None
        
        self.load_profile(with_stats=with_stats)

        self.sample_id = sample_id
        self.get_mask()

        self.calculate_rates(genome, cov_corr_without_gu=cov_corr_without_gu)

        if norm_without_gu:
            self.fac_norm = normalize(self.rat_p + self.rat_m, [self.mask_p[i] & genome.mask_seq[i] for i in range(len(genome.seq))] + [self.mask_m[i] & ~genome.mask_seq[i] for i in range(len(genome.seq))])
        else:
            self.fac_norm = normalize(self.rat_p + self.rat_m, [self.mask_p[i] for i in range(len(genome.seq))] + [self.mask_m[i] for i in range(len(genome.seq))])


    def load_profile(self, with_stats=False):
        with open(self.pkl_loc, 'rb') as f:
            if with_stats:
                self.mut_p, self.cov_p, self.mut_m, self.cov_m, self.stats = pickle.load(f)
            else:
                self.mut_p, self.cov_p, self.mut_m, self.cov_m, _ = pickle.load(f)
        
    def get_mask(self):
        self.mask_p = [(self.cov_p[i] >= self.min_cov) & (self.mut_p[i] >= self.min_mut) for i in range(len(self.cov_p))]
        self.mask_m = [(self.cov_m[i] >= self.min_cov) & (self.mut_m[i] >= self.min_mut) for i in range(len(self.cov_m))]

    def calculate_rates(self, genome, cov_corr_without_gu=True):
        # calculate parameters to correct for coverage bias
        if not cov_corr_without_gu:
            alpha, beta = normalize_cov(np.hstack(self.mut_p + self.mut_m)[np.hstack(self.mask_p + self.mask_m)], np.hstack(self.cov_p + self.cov_m)[np.hstack(self.mask_p + self.mask_m)], return_params=True)
        else:
            alpha, beta = normalize_cov(np.hstack(self.mut_p + self.mut_m)[np.hstack(self.mask_p + self.mask_m) & np.hstack(genome.mask_seq + [~i for i in genome.mask_seq])], np.hstack(self.cov_p + self.cov_m)[np.hstack(self.mask_p + self.mask_m) & np.hstack(genome.mask_seq + [~i for i in genome.mask_seq])], return_params=True)

        # apply
        self.rat_p = [(self.mut_p[i] + alpha) / (self.cov_p[i] + beta + alpha) for i in range(len(self.mut_p))]
        self.rat_m = [(self.mut_m[i] + alpha) / (self.cov_m[i] + beta + alpha) for i in range(len(self.mut_m))]     

    def plot_profile(self, chrom, strand, start, stop, genome, ax=None, plot_raw=False, cmap_loc='cmap.txt', **kwargs):

        cmap = np.array(pd.read_csv(cmap_loc, header=None))/255
        cmap = mplc.ListedColormap(cmap)
        norm = mplc.Normalize(vmin=0, vmax=1)
        m = cm.ScalarMappable(norm=norm, cmap=cmap)

        x = np.arange(start, stop)

        if plot_raw:
            f = 1
        else:
            f = self.fac_norm
        
        if strand == '+':
            y = self.rat_p[chrom][start:stop]
            clrs = m.to_rgba(self.fac_norm * y)
            clrs[~genome.mask_seq[chrom][start:stop]] = [230/255, 230/255, 230/255, 1]
            clrs[~self.mask_p[chrom][start:stop]] = [1, 0, 0, 1]
            if self.std_p is not None:
                error = self.std_p[chrom][start:stop]
        elif strand == '-':
            y = self.rat_m[chrom][start:stop]
            clrs = m.to_rgba(self.fac_norm * y)
            clrs[genome.mask_seq[chrom][start:stop]] = [230/255, 230/255, 230/255, 1] 
            clrs[~self.mask_m[chrom][start:stop]] = [1, 0, 0, 1]
            if self.std_m is not None:
                error = self.std_m[chrom][start:stop]
        
        if ax is not None:
            ax.bar(x, f*y, 1, color=clrs, **kwargs)
            if self.std_p is not None:
                ax.errorbar(x, f*y, yerr=f*error, color='k', ls='none')
        else:
            plt.bar(x, f*y, 1, color=clrs, **kwargs)
            if self.std_p is not None:
                plt.errorbar(x, f*y, yerr=f*error, color='k', ls='none')
            plt.show()


class PRO_Profile():
    """
    Class for handling PRO-seq (Precision Run-On Sequencing) profiling data.

    Parameters:
        - pkl_loc (str): Path to the pickle file containing PRO-seq profile data.
        - sample_id (str): Identifier for the PRO-seq profile sample.
        - genome (Genome): Genome object representing the genomic sequence.
        - norm (bool): Option to normalize PRO-seq profile depth (default: True).
        - with_stats (bool): Option to load statistics during initialization.

    Methods:
        - __init__: Initialize a PRO_Profile object with specified parameters.
        - load_profile: Load PRO-seq profile data from a pickle file.
        - normalize_depth: Normalize PRO-seq profile depth.
        - normalize_nt_content: Normalize PRO-seq profile nucleotide content.
        - plot_profile: Plot PRO-seq profile for a specified genomic region.
        - get_meta_profile_binned: Get binned meta-profile based on genomic annotation.
        - get_meta_profile: Get meta-profile based on genomic annotation.

    Attributes:
        - Various attributes for storing profile data and parameters.
    """
    
    def __init__(self, pkl_loc, sample_id, genome, norm=True, with_stats=False):
        self.pkl_loc = pkl_loc
        
        self.load_profile(with_stats=with_stats)

        if norm:
            self.normalize_depth()
            # self.normalize_nt_content(genome)

        self.sample_id = sample_id
        
    def load_profile(self, with_stats=False):
        with open(self.pkl_loc, 'rb') as f:
            if with_stats:
                self.pro_p, self.pro_m, self.stats = pickle.load(f)
            else:
                self.pro_p, self.pro_m, _ = pickle.load(f)

    def normalize_depth(self):
        self.cov_pm = np.sum(np.hstack(self.pro_p + self.pro_m)) / 1_000_000 # 3'ends per million
        self.pro_p = [i/self.cov_pm for i in self.pro_p]
        self.pro_m = [i/self.cov_pm for i in self.pro_m]

    def normalize_nt_content(self, genome):
        # problem for now is that at the beginning of a transcript this won't be accurate because it considers the pileup from the previous gene
        # but this only affects the very first nucleotide of every transcript, so it's probably okay just need to keep in mind
        
        self.pro_p = [self.pro_p[i]/genome.norm_p[i] for i in range(len(self.pro_p))]
        self.pro_m = [self.pro_m[i]/genome.norm_m[i] for i in range(len(self.pro_m))]

    def plot_profile(self, chrom, start, stop, genome, ax, avg_win=None, nts_runon=['A', 'C', 'G', 'U', 'T'], **kwargs):
        
        seq = genome.seq[chrom][start:stop]
        if avg_win is not None:
            if avg_win%2 != 0: # sliding window must be even number
                avg_win += 1
            x = np.arange(start, stop)[int(avg_win/2):-int(avg_win/2)+1]
            y_p =  moving_average(self.pro_p[chrom][start:stop], avg_win)
            y_m = -moving_average(self.pro_m[chrom][start:stop], avg_win)
            # colors
            clrs_p = [[230/255, 230/255, 230/255, 1] if i.upper() not in nts_runon else [240/255, 113/255, 21/255, 1] for i in seq[int(avg_win/2):-int(avg_win/2)+1]]
            clrs_m = [[230/255, 230/255, 230/255, 1] if i.upper() not in [reverse_complement(i) for i in nts_runon] else [158/255, 135/255, 94/255, 1] for i in seq[int(avg_win/2):-int(avg_win/2)+1]]
        else:
            x = np.arange(start, stop)
            y_p =  self.pro_p[chrom][start:stop]
            y_m = -self.pro_m[chrom][start:stop]
            # colors
            clrs_p = [[230/255, 230/255, 230/255, 1] if i.upper() not in nts_runon else [240/255, 113/255, 21/255, 1] for i in seq]
            clrs_m = [[230/255, 230/255, 230/255, 1] if i.upper() not in [reverse_complement(i) for i in nts_runon] else [158/255, 135/255, 94/255, 1] for i in seq]
        # plotting: y axis is nt-content corrected counts per million
        if ax is not None:
            ax.bar(x, y_p, 1, color=clrs_p, **kwargs)
            ax.bar(x, y_m, 1, color=clrs_m, **kwargs)
        else:
            plt.bar(x, y_p, 1, color=clrs_p, **kwargs)
            plt.bar(x, y_m, 1, color=clrs_m, **kwargs)
            plt.show()

    def get_meta_profile_binned(self, annotation, n_bins=100, surr=0.1):

        meta_binned = []
        chrom_dict = {'chrI':0, 'chrII':1, 'chrIII':2, 'chrIV':3, 'chrV':4, 'chrVI':5, 'chrVII':6, 'chrVIII':7, 'chrIX':8, 'chrX':9, 'chrXI':10, 'chrXII':11, 'chrXIII':12, 'chrXIV':13, 'chrXV':14, 'chrXVI':15, 'chrmt':16, 'KanMX': 17}
        
        for gene in annotation.intron_less + annotation.intron_cont:
            feat = annotation.db[gene]
            start = [feat.start]
            end = [feat.end]
                    
            extend = int(np.round(surr*(end[0]-start[0]+1)))

            # skip edge cases
            if start[0]-extend-1 < 0 or end[0]+extend >= len(self.pro_p[chrom_dict[feat.seqid]]):
                continue

            if feat.strand == '+':
                signal_y = self.pro_p[chrom_dict[feat.seqid]][start[0]-extend-1:end[0]+extend]
            elif feat.strand == '-':
                signal_y = (self.pro_m[chrom_dict[feat.seqid]][start[0]-extend-1:end[0]+extend])[::-1]

            norm_x, norm_y, bin_size = binned_signal(signal_y, n_bins)
            meta_binned.append(norm_y/np.sum(norm_y))


        # return([i for i in meta_binned if i == i], norm_x, bin_size)
        return(np.array(meta_binned), norm_x, bin_size)

    def get_meta_profile(self, annotation, surr=(50, 100)):

        meta_start = []
        meta_stop  = []
        chrom_dict = {'chrI':0, 'chrII':1, 'chrIII':2, 'chrIV':3, 'chrV':4, 'chrVI':5, 'chrVII':6, 'chrVIII':7, 'chrIX':8, 'chrX':9, 'chrXI':10, 'chrXII':11, 'chrXIII':12, 'chrXIV':13, 'chrXV':14, 'chrXVI':15, 'chrmt':16, 'KanMX': 17}
        
        for gene in annotation.intron_less + annotation.intron_cont:
            feat = annotation.db[gene]
            start = [feat.start]
            end = [feat.end]
                    
            # skip edge cases
            if start[0]-max(surr)-1 < 0 or end[0]+max(surr) >= len(self.pro_p[chrom_dict[feat.seqid]]):
                continue

            if feat.strand == '+':
                signal_start = self.pro_p[chrom_dict[feat.seqid]][start[0]-surr[0]:start[0]+surr[1]]
                signal_stop  = self.pro_p[chrom_dict[feat.seqid]][end[0]-surr[1]:end[0]+surr[0]]
            elif feat.strand == '-':
                signal_stop  = (self.pro_p[chrom_dict[feat.seqid]][start[0]-surr[0]:start[0]+surr[1]])[::-1]
                signal_start = (self.pro_p[chrom_dict[feat.seqid]][end[0]-surr[1]:end[0]+surr[0]])[::-1]

            meta_start.append(signal_start/np.sum(signal_start + signal_stop))
            meta_stop.append(signal_stop/np.sum(signal_stop + signal_start))


        return(meta_start, meta_stop)
        

class Genome():
    """
    Class for handling genomic sequence data.

    Parameters:
        - filename (str): Path to the fasta file containing genomic sequence data.
        - coords (tuple): Genomic coordinates for subsetting sequence data (default: None).
        - reverse (bool): Option to reverse the genomic sequence, only applies if subsetting (default: False).

    Methods:
        - __init__: Initialize a Genome object with genomic sequence data.
        - mask_A: Return a mask for 'A' nucleotides.
        - mask_G: Return a mask for 'G' nucleotides.
        - mask_C: Return a mask for 'C' nucleotides.
        - mask_U: Return a mask for 'U' nucleotides.

    Attributes:
        - Various attributes for storing genomic sequence data and parameters.
    """

    def __init__(self, filename, coords=None, reverse=False):
        self.seq = []
        self.reverse=reverse
        self.coords=coords
        for i in SeqIO.parse(filename, 'fasta'):
            self.seq.append(str(i.seq).upper())

        if coords is not None:
            self.seq = self.seq[coords[0]][coords[1]:coords[2]]
            if not reverse:
                self.mask_seq = np.array([True if i in ['A', 'C'] else False for i in self.seq])
            else:
                self.mask_seq = np.array([True if i in ['A', 'C'] else False for i in reverse_complement(self.seq)])
                self.seq = reverse_complement(self.seq)
        
        else:        
            self.mask_seq = [np.array([True if i in ['A', 'C'] else False for i in chrom]) for chrom in self.seq]

            # calculations for PROseq normalization
            self.norm_p = [c_freq(i) for i in self.seq]
            self.norm_m = [c_freq(reverse_complement(i))[::-1] for i in self.seq]    

    def mask_A(self):
        if self.coords is None:
            return([np.array([True if i.upper() == 'A' else False for i in chrom]) for chrom in self.seq])
        else:
            return(np.array([True if i.upper() == 'A' else False for i in self.seq]))
    def mask_G(self):
        if self.coords is None:
            return([np.array([True if i.upper() == 'G' else False for i in chrom]) for chrom in self.seq])
        else:
            return(np.array([True if i.upper() == 'G' else False for i in self.seq]))
    def mask_C(self):
        if self.coords is None:
            return([np.array([True if i.upper() == 'C' else False for i in chrom]) for chrom in self.seq])
        else:
            return(np.array([True if i.upper() == 'C' else False for i in self.seq]))
    def mask_U(self):
        if self.coords is None:       
            return([np.array([True if i.upper() in ['U', 'T'] else False for i in chrom]) for chrom in self.seq])
        else:
            return(np.array([True if i.upper() in ['U', 'T'] else False for i in self.seq]))

class Annotation():
    """
    Class for handling genomic annotation data.

    Parameters:
        - ann_loc (str): Path to the .gff file containing genomic annotation data.

    Methods:
        - __init__: Initialize an Annotation object with genomic annotation data.

    Attributes:
        - Various attributes for storing genomic annotation data.
    """

    def __init__(self, ann_loc):
        self.db = gffutils.create_db(ann_loc, dbfn='ann.db', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)

        # select certain transcripts
        # intron containing and intronless mRNAs (exclude weird introns)
        intron_less = []
        intron_cont = []
        for i in self.db.features_of_type('mRNA'):
            p = [child.featuretype for child in self.db.children(i.id)]
            if 'intron' in p and not any(x in p for x in ['five_prime_UTR_intron', 'telomeric_repeat', 'intein_encoding_region', 'plus_1_translational_frameshift', 'uORF']):
                if 'conditions' in i.attributes:
                    if 'YPD' in i.attributes['conditions']:
                        intron_cont.append(i.id)
                else:
                    intron_cont.append(i.id)
            elif 'intron' not in p and not any(x in p for x in ['five_prime_UTR_intron', 'telomeric_repeat', 'intein_encoding_region', 'plus_1_translational_frameshift', 'uORF']):
                if 'conditions' in i.attributes:
                    if 'YPD' in i.attributes['conditions']:
                        intron_less.append(i.id)
                else:
                    intron_less.append(i.id)

        print(f"intron-containing transcripts: {len(intron_cont)}")
        print(f"intron-less transcripts: {len(intron_less)}")

        # structured RNAs (exclude intron-containing snRNAs)
        structured = []
        for i in self.db.features_of_type(('ncRNA_gene', 'snoRNA_gene', 'snRNA_gene', 'telomerase_RNA_gene')):
            p = [child.featuretype for child in self.db.children(i.id)]
            if 'intron' not in p:
                if 'conditions' in i.attributes:
                    if 'YPD' in i.attributes['conditions']:
                        structured.append(i.id)
                else:
                    structured.append(i.id)
        print(f"structured RNA transcripts: {len(structured)}")
        self.intron_cont = intron_cont
        self.intron_less = intron_less
        self.structured = structured

class HDP():
    """
    Class for handling HDProbe data.

    Parameters:
        - filename (str): Path to the file containing HDP data.
        - sample_id (str): Identifier for the HDP sample.
        - thresh_col (str): Threshold column for significance (default: 'padj_BF').

    Methods:
        - __init__: Initialize an HDP object with specified parameters.
        - get_distance_to_nearest_hit: Calculate distances to nearest significant hits.
        - annotate_hits: Annotate hits based on genomic annotation.
        - erupt: Volcano plot for significant hits.
        - get_hdp_meta: Extract meta-profiles from significant hits.

    Attributes:
        - Various attributes for storing HDP data and parameters.
    """

    def __init__(self, filename, sample_id, thresh_col='padj_BF'):
        self.hdp = pd.read_csv(filename)
        self.sample_id = sample_id
        # filter to save memory
        self.hdp = self.hdp[self.hdp[thresh_col] <= 0.5]
        self.thresh_col = thresh_col

        self.get_distance_to_nearest_hit(thresh_col=thresh_col)

    def get_distance_to_nearest_hit(self, sig_threshold=0.05, thresh_col='padj_BF'):
        dists = []
        hdp_f = self.hdp[self.hdp[thresh_col] <= sig_threshold].copy()
        for idx, hit in hdp_f.iterrows():
            all_dists = np.abs(hdp_f[np.logical_and(hdp_f['chromosome']==hit['chromosome'], hdp_f['strand']==hit['strand'])]['position'] - hit['position'])
            if len(all_dists) == 1:
                dists.append(999999999)
            else:
                dists.append(np.partition(all_dists, 1)[1])

        hdp_f.loc[:, 'min_dist'] = dists
        self.hdp = hdp_f

    def annotate_hits(self, annotation, surr=50):
        ann = []
        chrom_dict = {0:'chrI', 1:'chrII', 2:'chrIII', 3:'chrIV', 4:'chrV', 5:'chrVI', 6:'chrVII', 7:'chrVIII', 8:'chrIX', 9:'chrX', 10:'chrXI', 11:'chrXII', 12:'chrXIII', 13:'chrXIV', 14:'chrXV', 15:'chrXVI', 16:'chrmt', 17:'KanMX'}
        for idx, hit in self.hdp.iterrows():
            features = feature_search(annotation.db, chrom_dict[hit['chromosome']-1], hit['position'], surr=surr, strand=hit['strand'])
            # if there is more than one feature mapping to the position: 
            if len(features) > 1:
                features = [i for i in features if 'orf_classification=Dubious' not in str(i).split(';')]
            if not features:
                features = [{'ID':['???']}]

            ann.append(features[0]['ID'][0])

        self.hdp['annotation'] = ann

    def erupt(self, axs, cmap_loc='cmap.txt', **kwargs):

        cmap = np.flipud(np.array(pd.read_csv(cmap_loc, header=None))/255)
        cmap = mplc.ListedColormap(cmap)
        norm = mplc.Normalize(vmin=0, vmax=3)
        m = cm.ScalarMappable(norm=norm, cmap=cmap)

        clrs = m.to_rgba(np.log10(self.hdp['min_dist']))

        axs.scatter(self.hdp['difference'], -np.log10(self.hdp[self.thresh_col]), 10, clrs, **kwargs)

    
    def get_hdp_meta(self, annotation, cov_p, cov_m, thresh_col='padj_BH', sig_threshold=0.01, diff_threshold=0.8, binsize=0.02, n_bins_cov=100):
        # where in the gene is the hit?
        abso_all = []
        anno_all = []
        rela_all = []
        rcov_all = []
        chrom_dict = {'chrI':0, 'chrII':1, 'chrIII':2, 'chrIV':3, 'chrV':4, 'chrVI':5, 'chrVII':6, 'chrVIII':7, 'chrIX':8, 'chrX':9, 'chrXI':10, 'chrXII':11, 'chrXIII':12, 'chrXIV':13, 'chrXV':14, 'chrXVI':15, 'chrmt':16, 'KanMX': 17}
        for idx, row in self.hdp.iterrows():
            if row[thresh_col] <= sig_threshold and row['annotation'] != '???' and abs(row['difference']) >= diff_threshold:
                # determine absolute position of hit in transcript
                leng = abs(annotation.db[row['annotation']].start - annotation.db[row['annotation']].end)
                if row['strand'] == '+':
                    abso = row['position'] - annotation.db[row['annotation']].start
                    cov  = cov_p[chrom_dict[annotation.db[row['annotation']].seqid]][annotation.db[row['annotation']].start-int(0.5*leng):annotation.db[row['annotation']].end+int(0.5*leng)] 
                elif row['strand'] == '-':
                    abso = annotation.db[row['annotation']].end - row['position']
                    cov  = cov_m[chrom_dict[annotation.db[row['annotation']].seqid]][annotation.db[row['annotation']].start-int(0.5*leng):annotation.db[row['annotation']].end+int(0.5*leng)]
                    cov = cov[::-1]

                # determine relative position of hit in transcript
                rela = abso/leng
                _, covy, _ = binned_signal(cov, n_bins_cov)

                abso_all.append(abso)
                anno_all.append(row['annotation'])
                rela_all.append(rela)
                rcov_all.append(covy)

        bins_hit = np.arange(-0.5, 1.5+binsize, binsize)
        bins_cov, binsize_cov = np.linspace(-0.5, 1.5+binsize, num=n_bins_cov, retstep=True)
        hist_hit, _ = np.histogram(rela_all, bins=bins_hit, density=False)
        hist_cov = np.median(rcov_all, axis=0)

        return(bins_hit[:-1]+0.5*binsize, hist_hit, bins_cov, hist_cov, rela_all, anno_all)


class Targeted_DMS_Profile():
    """
    Class for handling targeted DMS profiling data.

    Parameters:
        - pkl_loc (str): Path to the pickle file containing targeted DMS profile data.
        - sample_id (str): Identifier for the targeted DMS profile sample.
        - genome (Genome): Genome object representing the genomic sequence.
        - min_cov (int): Minimum coverage threshold for filtering data (default: 1000).
        - min_mut (int): Minimum mutation threshold for filtering data (default: 1).
        - cov_corr_individual (bool): Option for individual coverage correction.
        - reverse (bool): Option to reverse the genomic sequence (default: False).
        - use_G (bool): Option to include 'G' bases in analysis.
        - use_U (bool): Option to include 'U' bases in analysis (default: True).
        - per_exclude (int): Upper percentile of positions to exclude from normalization (default: 10).
        - with_stats (bool): Option to load statistics during initialization.

    Methods:
        - __init__: Initialize a Targeted_DMS_Profile object with specified parameters.
        - load_profile: Load targeted DMS profile data from a pickle file.
        - get_mask: Create masks based on coverage and mutation thresholds.
        - calculate_rates: Calculate mutation rates with optional coverage correction.
        - plot_profile: Plot targeted DMS profile for a specified genomic region.

    Attributes:
        - Various attributes for storing profile data and parameters.
    """

    def __init__(self, pkl_loc, sample_id, genome, min_cov=1000, min_mut=1, cov_corr_individual=True, reverse=False, use_G=False, use_U=True, per_exclude=10, with_stats=False):
        # parameters
        self.min_cov = min_cov
        self.min_mut = min_mut

        self.pkl_loc = pkl_loc
        
        self.load_profile(reverse=reverse, with_stats=with_stats)

        self.sample_id = sample_id
        self.get_mask()
        self.calculate_rates(genome, cov_corr_individual=cov_corr_individual)

        self.fac_norm = {}

        self.fac_norm['A'] = normalize([self.rat], [genome.mask_A()], per_exclude=per_exclude)
        self.fac_norm['C'] = normalize([self.rat], [genome.mask_C()], per_exclude=per_exclude)
        self.fac_norm['G'] = normalize([self.rat], [genome.mask_G()], per_exclude=per_exclude) if use_G else 1
        self.fac_norm['U'] = normalize([self.rat], [genome.mask_U()], per_exclude=per_exclude) if use_U else 1
        self.fac_norm['T'] = self.fac_norm['U']
        self.vec_norm = np.array([self.fac_norm[nt.upper()] for nt in genome.seq])

    def load_profile(self, reverse=False, with_stats=False):
        with open(self.pkl_loc, 'rb') as f:
            if with_stats:
                self.mut, self.cov, self.stats = pickle.load(f)
            else:
                self.mut, self.cov, _ = pickle.load(f)
        if reverse:
            self.mut = self.mut[::-1]
            self.cov = self.cov[::-1]

        
    def get_mask(self):
        self.mask = (self.cov >= self.min_cov) & (self.mut >= self.min_mut)

    def calculate_rates(self, genome, cov_corr_individual=True):
        # calculate parameters to correct for coverage bias, each nt separately
        if cov_corr_individual:
            alpha_a, beta_a = normalize_cov(self.mut[self.mask & genome.mask_A()], self.cov[self.mask & genome.mask_A()], return_params=True)
            alpha_c, beta_c = normalize_cov(self.mut[self.mask & genome.mask_C()], self.cov[self.mask & genome.mask_C()], return_params=True)
            alpha_u, beta_u = normalize_cov(self.mut[self.mask & genome.mask_U()], self.cov[self.mask & genome.mask_U()], return_params=True)
            alpha_g, beta_g = normalize_cov(self.mut[self.mask & genome.mask_G()], self.cov[self.mask & genome.mask_G()], return_params=True)

            # apply
            corr_vec_alpha = []
            corr_vec_beta = []
            for nt in genome.seq:
                if nt.upper() == 'A':
                     corr_vec_alpha.append(alpha_a)
                     corr_vec_beta.append(beta_a)
                if nt.upper() == 'C':
                     corr_vec_alpha.append(alpha_c)
                     corr_vec_beta.append(beta_c)
                if nt.upper() == 'G':
                     corr_vec_alpha.append(alpha_g)
                     corr_vec_beta.append(beta_g)
                if nt.upper() in ['U', 'T']:
                     corr_vec_alpha.append(alpha_u)
                     corr_vec_beta.append(beta_u)

            self.rat = (self.mut + np.array(corr_vec_alpha)) / (self.cov + np.array(corr_vec_beta))

        else:
            alpha, beta = normalize_cov(self.mut[self.mask & (genome.mask_A() | genome.mask_C() | genome.mask_U())], self.cov[self.mask & (genome.mask_A() | genome.mask_C() | genome.mask_U())], return_params=True)
            self.rat = (self.mut + alpha) / (self.cov + beta)

    def plot_profile(self, start, stop, genome, ax, raw_rats=False, cmap_loc='cmap.txt', **kwargs):

        cmap = np.array(pd.read_csv(cmap_loc, header=None))/255
        cmap = mplc.ListedColormap(cmap)
        norm = mplc.Normalize(vmin=0, vmax=1)
        m = cm.ScalarMappable(norm=norm, cmap=cmap)

        x = np.arange(start, stop)
        y = self.rat[start:stop]
        
        clrs = m.to_rgba(y * self.vec_norm[start:stop])

        if self.fac_norm['G'] == 1: 
            clrs[genome.mask_G()[start:stop] | ~self.mask[start:stop]] = [230/255, 230/255, 230/255, 1]
        else:
            clrs[~self.mask[start:stop]] = [230/255, 230/255, 230/255, 1]

        if hasattr(self, 'std'):
            error = self.std[start:stop]  
        
        if not raw_rats:
            ax.bar(x, y * self.vec_norm[start:stop], 1, color=clrs, **kwargs)
            if hasattr(self, 'std'):
                ax.errorbar(x, y * self.vec_norm[start:stop], yerr=error * self.vec_norm[start:stop], color='k', ls='none')
        else:
            ax.bar(x, y, 1, color=clrs, **kwargs)
            if hasattr(self, 'std'):
                ax.errorbar(x, y, yerr=error, color='k', ls='none')

        return(x, y, self.vec_norm[start:stop])


class DMS_Profile_just_rea():

    """
    Class for handling DMS profiling data (only reactivity information).

    Parameters:
        - rat (str): Reactivity data for DMS profiling.
        - sample_id (str): Identifier for the DMS profile sample.
        - use_G (bool): Option to include 'G' bases in analysis.
        - use_U (bool): Option to include 'U' bases in analysis (default: True).
        - per_exclude (int): Upper percentile of positions to exclude from normalization (default: 10).

    Methods:
        - __init__: Initialize a DMS_Profile_just_rea object with specified parameters.
        - plot_profile: Plot DMS profile for a specified genomic region.

    Attributes:
        - Various attributes for storing reactivity data and parameters.
    """

    def __init__(self, rat, sample_id, use_G=False, use_U=True, per_exclude=10):
        # parameters
        self.sample_id = sample_id
        self.rat = rat

    def plot_profile(self, start, stop, genome, ax, raw_rats=False, cmap_loc='cmap.txt', **kwargs):

        cmap = np.array(pd.read_csv(cmap_loc, header=None))/255
        cmap = mplc.ListedColormap(cmap)
        norm = mplc.Normalize(vmin=0, vmax=0.5)
        m = cm.ScalarMappable(norm=norm, cmap=cmap)

        x = np.arange(start, stop)
        y = self.rat[start:stop]
        
        clrs = m.to_rgba(y)

        clrs[genome.mask_G()[start:stop]] = [230/255, 230/255, 230/255, 1]
       
        ax.bar(x, y, 1, color=clrs, **kwargs)

        return(x, y)


def combine_profiles(profiles, sample_id, genome, norm_without_gu=False, cov_corr_without_gu=True, min_cov=700):
    """
    Combine multiple DMS profiles into a single profile.

    Parameters:
    - profiles (List[DMS_Profile]): List of DMS_Profile objects to be combined.
    - sample_id (str): Sample ID for the new combined profile.
    - genome (Genome): Genome object.
    - norm_without_gu (bool): Flag to indicate whether to normalize without considering G and U nucleotides (default is False).
    - cov_corr_without_gu (bool): Flag to indicate whether to perform coverage correction without considering G and U nucleotides (default is True).
    - min_cov (int): Minimum coverage threshold (default is 700).

    Returns:
    - DMS_Profile: Combined DMS profile.
    """

    assert isinstance(profiles, list), "Input must be list of DMS_Profile objects."

    new = deepcopy(profiles[0])
    delattr(new, 'pkl_loc')
    new.sample_id = sample_id
    new.min_cov = min_cov

    new.cov_p = [np.sum([i.cov_p[chrom] for i in profiles], axis=0) for chrom in range(len(profiles[0].cov_p))]
    new.cov_m = [np.sum([i.cov_m[chrom] for i in profiles], axis=0) for chrom in range(len(profiles[0].cov_m))]
    new.mut_p = [np.sum([i.mut_p[chrom] for i in profiles], axis=0) for chrom in range(len(profiles[0].mut_p))]
    new.mut_m = [np.sum([i.mut_m[chrom] for i in profiles], axis=0) for chrom in range(len(profiles[0].mut_m))]

    new.get_mask()
    new.calculate_rates(genome, cov_corr_without_gu=cov_corr_without_gu)

    # calculate error based on replicates
    new.std_p = []
    new.std_m = []
    for chrom in range(len(genome.seq)):
        rat_ps = []
        rat_ms = []
        for profile in profiles:
            rat_p_corr = profile.rat_p[chrom].copy()
            rat_p_corr[~profile.mask_p[chrom]] = np.nan
            rat_m_corr = profile.rat_m[chrom].copy()
            rat_m_corr[~profile.mask_m[chrom]] = np.nan
            rat_ps.append(rat_p_corr)
            rat_ms.append(rat_m_corr)
        new.std_p.append(np.nanstd(rat_ps, axis=0))
        new.std_m.append(np.nanstd(rat_ms, axis=0))

    if norm_without_gu:
        new.fac_norm = normalize(new.rat_p + new.rat_m, [new.mask_p[i] & genome.mask_seq[i] for i in range(len(genome.seq))] + [new.mask_m[i] & ~genome.mask_seq[i] for i in range(len(genome.seq))])
    else:
        new.fac_norm = normalize(new.rat_p + new.rat_m, [new.mask_p[i] for i in range(len(genome.seq))] + [new.mask_m[i] for i in range(len(genome.seq))])

    return(new)


def targeted_combine_profiles(profiles, sample_id, genome, cov_corr_individual=True, per_exclude=10, use_G=False, use_U=True, min_cov=700):
    """
    Combine targeted DMS profiles into a single profile.

    Parameters:
    - profiles (List[Targeted_DMS_Profile]): List of Targeted_DMS_Profile objects to be combined.
    - sample_id (str): Sample ID for the new combined profile.
    - genome (Genome): Genome object.
    - cov_corr_individual (bool): Flag to indicate whether to perform coverage correction individually (default is True).
    - per_exclude (int): Percentage of values to exclude during normalization (default is 10).
    - use_G (bool): Flag to include G nucleotides in normalization (default is False).
    - use_U (bool): Flag to include U nucleotides in normalization (default is True).
    - min_cov (int): Minimum coverage threshold (default is 700).

    Returns:
    - DMS_Profile: Combined targeted DMS profile.
    """

    assert isinstance(profiles, list), "Input must be list of DMS_Profile objects."

    new = deepcopy(profiles[0])
    delattr(new, 'pkl_loc')
    new.sample_id = sample_id
    new.min_cov = min_cov

    new.cov = np.sum([i.cov for i in profiles], axis=0)
    new.mut = np.sum([i.mut for i in profiles], axis=0)

    new.get_mask()
    new.calculate_rates(genome, cov_corr_individual=cov_corr_individual)
    
    # calculate error based on replicates, only include nucleotide if it reaches min_cov in individual replicate
    rats = []
    for profile in profiles:
        rat_corr = profile.rat.copy()
        rat_corr[~profile.mask] = np.nan
        rats.append(rat_corr)
    new.std = np.nanstd(rats, axis=0)

    new.fac_norm = {}
    new.fac_norm['A'] = normalize(new.rat, genome.mask_A(), per_exclude=per_exclude)
    new.fac_norm['C'] = normalize(new.rat, genome.mask_C(), per_exclude=per_exclude)
    new.fac_norm['G'] = normalize(new.rat, genome.mask_G(), per_exclude=per_exclude) if use_G else 1
    new.fac_norm['U'] = normalize(new.rat, genome.mask_U(), per_exclude=per_exclude) if use_U else 1
    new.fac_norm['T'] = new.fac_norm['U']
    new.vec_norm = np.array([new.fac_norm[nt.upper()] for nt in genome.seq])

    return(new)

def pro_combine_profiles(profiles, sample_id, genome, norm=True):
    """
    Combine profiles for PROseq data.

    Parameters:
    - profiles (List[PRO_Profile]): List of PRO_Profile objects.
    - sample_id (str): Sample ID for the new combined profile.
    - genome (Genome): Genome object.
    - norm (bool): Flag to indicate whether to perform normalization (default is True).

    Returns:
    - PRO_Profile: Combined PROseq profile.
    """

    assert isinstance(profiles, list), "Input must be list of DMS_Profile objects."

    new = deepcopy(profiles[0])
    delattr(new, 'pkl_loc')
    new.sample_id = sample_id

    new.pro_p = [np.sum([i.pro_p[chrom] for i in profiles], axis=0) for chrom in range(len(profiles[0].pro_p))]
    new.pro_m = [np.sum([i.pro_m[chrom] for i in profiles], axis=0) for chrom in range(len(profiles[0].pro_m))]

    if norm:
        new.normalize_depth()
        new.normalize_nt_content(genome)

    return(new)

def normalize(arr, mask=None, per_exclude=10):
    """
    Normalize an array of reactivity values.

    Parameters:
    - arr (numpy.ndarray): Array of values to be normalized.
    - mask (numpy.ndarray): Mask for values to be considered in normalization (default is None).
    - per_exclude (int): Upper percentile of values to exclude during normalization (default is 10).

    Returns:
    - float: Normalization factor.
    """

    # sort and remove nan
    arr_sort = []

    if mask is not None:
        for msk, chrom in zip(mask, arr):
            arr_sort += list(chrom[msk])
    else:
        for chrom in arr:
            arr_sort += list(chrom)

    arr_sort = np.array(arr_sort)
    arr_sort = np.sort(arr_sort[arr_sort==arr_sort])

    if len(arr_sort) == 0:
        return(1)

    # define values excluded for normalization
    IQR = iqr(arr_sort)
    PER = arr_sort[arr_sort.shape[0] - 1 - (len(arr_sort)//per_exclude if len(arr_sort) > 100 else len(arr_sort)//5)]
    threshold = np.max([1.5*IQR, PER])
    good = (arr_sort < threshold)
    
    # calculate normalization factor (mean of top 10 percent -> length according to ORIGINAL vector)
    # I think length should be according to the new vector (after excluding values), so potentially change this in the future
    top_rea = arr_sort[good][(-len(arr_sort)//10)+1:]
    fac_norm = 1/np.mean(top_rea)
    
    return(fac_norm)



def write_HDProbe(out_name, profiles, seq, min_cov=500, max_rate=0.3):
    """
    Write HDProbe format file.

    Parameters:
    - out_name (str): Output file name.
    - profiles (List[DMS_Profile]): List of DMS_Profile objects.
    - seq (str): Genomic sequence.
    - min_cov (int): Minimum coverage threshold (default is 500).
    - max_rate (float): Maximum mutation rate (default is 0.3).

    Returns:
    - None
    """
    
    assert isinstance(profiles, list), "Input must be list of DMS_Profile objects."
    
    with open(out_name, 'w') as f:
        # header
        f.write('\t'.join(['sample_name', 'chromosome', 'genomic_position', 'nt', 'strand', 'n_mutations', 'n_reads' + '\n']))
        
        for profile in tqdm(profiles):    
            for chrom in range(len(seq)):
                for nt in range(len(seq[chrom])):
                    # filter for coverage and SNPs
                    if profile.cov_p[chrom][nt] >= min_cov and profile.mut_p[chrom][nt]/profile.cov_p[chrom][nt] <= max_rate:
                        f.write('\t'.join([profile.sample_id, str(chrom+1), str(nt+1), seq[chrom][nt], '+', str(profile.mut_p[chrom][nt]), str(profile.cov_p[chrom][nt]) + '\n']))
                    if profile.cov_m[chrom][nt] >= min_cov and profile.mut_m[chrom][nt]/profile.cov_m[chrom][nt] <= max_rate:
                        f.write('\t'.join([profile.sample_id, str(chrom+1), str(nt+1), reverse_complement(seq[chrom][nt]), '-', str(profile.mut_m[chrom][nt]), str(profile.cov_m[chrom][nt]) + '\n']))


def targeted_write_HDProbe(out_name, profiles, seq, min_cov=500):
    """
    Write HDProbe format file for targeted profiles.

    Parameters:
    - out_name (str): Output file name.
    - profiles (List[Targeted_DMS_Profile]): List of Targeted_DMS_Profile objects.
    - seq (str): Target-specific sequence.
    - min_cov (int): Minimum coverage threshold (default is 500).

    Returns:
    - None
    """
    
    assert isinstance(profiles, list), "Input must be list of DMS_Profile objects."
    
    with open(out_name, 'w') as f:
        # header
        f.write('\t'.join(['sample_name', 'chromosome', 'genomic_position', 'nt', 'strand', 'n_mutations', 'n_reads' + '\n']))
        for profile in profiles:                            
            for nt in range(len(seq)):
                # filter for coverage and SNPs
                if profile.cov[nt] >= min_cov:
                    f.write('\t'.join([profile.sample_id, '11', str(nt+1), seq[nt], '+', str(profile.mut[nt]), str(profile.cov[nt]) + '\n']))
                
        

def gini_for_transcripts(profile1, profile2, transcript_list, genome, annotation, min_nt_gini=30, min_frac_of_transcript_gini=0.1):
    """
    Calculate Gini coefficients and correlations for transcripts.

    Parameters:
    - profile1 (DMS_Profile): First DMS_Profile object.
    - profile2 (DMS_Profile): Second DMS_Profile object.
    - transcript_list (List[str]): List of transcript IDs.
    - genome (Genome): Genome object.
    - annotation: Annotation object.
    - min_nt_gini (int): Minimum nucleotides for Gini calculation (default is 30).
    - min_frac_of_transcript_gini (float): Minimum fraction of transcript covered for Gini calculation (default is 0.1).

    Returns:
    - Tuple[Dict[str, Dict[str, Tuple[float, float]]], Dict[str, Dict[str, float]]]: Gini coefficients and correlations for transcripts.
    """

    chrom_dict = {'chrI':0, 'chrII':1, 'chrIII':2, 'chrIV':3, 'chrV':4, 'chrVI':5, 'chrVII':6, 'chrVIII':7, 'chrIX':8, 'chrX':9, 'chrXI':10, 'chrXII':11, 'chrXIII':12, 'chrXIV':13, 'chrXV':14, 'chrXVI':15, 'chrmt':16, 'KanMX': 17}
    
    gini_transcript = defaultdict(dict)
    corr_transcript = defaultdict(dict)
    
    for transcript in transcript_list:
        for child in annotation.db.children(transcript):
            if child.featuretype in ['CDS', 'noncoding_exon', 'intron']:
                if child.strand == '+':
                    rat1 = profile1.rat_p[chrom_dict[child.seqid]][child.start-1:child.end]
                    rat2 = profile2.rat_p[chrom_dict[child.seqid]][child.start-1:child.end]
                    msk = profile1.mask_p[chrom_dict[child.seqid]][child.start-1:child.end] & profile2.mask_p[chrom_dict[child.seqid]][child.start-1:child.end] & genome.mask_seq[chrom_dict[child.seqid]][child.start-1:child.end]
                elif child.strand == '-':
                    rat1 = profile1.rat_m[chrom_dict[child.seqid]][child.start-1:child.end]
                    rat2 = profile2.rat_m[chrom_dict[child.seqid]][child.start-1:child.end]
                    msk = profile1.mask_m[chrom_dict[child.seqid]][child.start-1:child.end] & profile2.mask_m[chrom_dict[child.seqid]][child.start-1:child.end] & ~genome.mask_seq[chrom_dict[child.seqid]][child.start-1:child.end]


                if (np.sum(msk) >= min_nt_gini) & (np.sum(msk)/(child.end-child.start+1) >= min_frac_of_transcript_gini):
                    gini_transcript[child.featuretype][transcript] = (gini(rat1[msk]), gini(rat2[msk]))
                    corr_transcript[child.featuretype][transcript] = pearsonr(rat1[msk], rat2[msk]).statistic

            else:
                print(child.featuretype)

    return(gini_transcript, corr_transcript)

def bulk_2_target(profile, chrom, strand, coords, genome):
    """
    Convert bulk DMS profile to targeted DMS profile.

    Parameters:
    - profile (DMS_Profile): Bulk DMS_Profile object.
    - chrom (int): Chromosome index.
    - strand (str): Strand information.
    - coords (Tuple[int, int]): Genomic coordinates.
    - genome (Genome): Genome object.

    Returns:
    - Targeted_DMS_Profile: Targeted DMS profile.
    """
    with open('temp_prof.pkl', 'wb') as f:
        if strand == '+':
            pickle.dump((profile.mut_p[chrom][coords[0]:coords[1]], profile.cov_p[chrom][coords[0]:coords[1]], 0), f)
            reverse = False
        else:
            pickle.dump((profile.mut_m[chrom][coords[0]:coords[1]], profile.cov_m[chrom][coords[0]:coords[1]], 0), f)
            reverse = True

    profile_t = Targeted_DMS_Profile('temp_prof.pkl', 'ddd', genome, reverse=reverse, min_cov=1000, min_mut=1, cov_corr_individual=True, use_G=False, use_U=False, per_exclude=5, with_stats=False)
    return(profile_t)


def feature_search(db, chrom, pos, surr=0, strand=None):
    return([i for i in db.region(seqid=chrom, start=pos-surr, end=pos+surr, strand=strand, 
                                     completely_within=False, featuretype=['gene', 'ncRNA_gene', 'pseudogene',
                                     'rRNA_gene', 'snRNA_gene', 'snoRNA_gene', 'tRNA_gene',
                                     'telomerase_RNA_gene', 'transposable_element_gene'])])

    """
    Search for genomic features around a specific position.

    Parameters:
    - db: Annotation database (gffutils object).
    - chrom (str): Chromosome or sequence identifier.
    - pos (int): Genomic position.
    - surr (int): Surrounding region size for the search (default is 0).
    - strand (str): Strand information (default is None).

    Returns:
    - List[Feature]: List of genomic features found in the specified region.
    """


















