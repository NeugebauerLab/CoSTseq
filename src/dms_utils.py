from src.bam_utils import read_pair_generator, reverse_complement, complement
import pysam
import numpy as np
from collections import Counter
from multiprocessing import Pool
from Bio import SeqIO, Seq
from scipy.stats import mannwhitneyu, pearsonr
import pandas as pd
import time
import re
import subprocess
import tempfile
import os.path
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors as mplc
import matplotlib.pyplot as plt

def coverage_mutations_PE(file, genome_file, n_cpu=4, clip_filter=9999999, clip_filter_3=True, clip_3=0, return_read_stats=False, allow_polyA_clipping=False):
    """
    Parallel wrapper for coverage_mutations_chrom_PE.

    Args:
        file (str): Path to the BAM file.
        genome_file (str): Path to the genome file in FASTA format.
        n_cpu (int): Number of CPU cores for parallel processing.
        clip_filter (int): Maximum allowed number of clipped bases.
        clip_filter_3 (bool): Whether to filter based on 3'-end clipping.
        clip_3 (int): Number of bases to remove from the 3'-end.
        return_read_stats (bool): Whether to return read statistics (generates large output file).
        allow_polyA_clipping (bool): Whether to allow poly(A) clipping.

    Returns:
        tuple: A tuple containing vectors with np.arrays, one for each chromosome for: mutation count plus strand, coverage plus strand, mutation count minus srand, coverage minus strand, stats (optional)
    """

    start = time.time()

    # load the genome    
    chroms = []
    for i in SeqIO.parse(genome_file, 'fasta'):
        chroms.append(i)
    
    vec_chroms = []
    
    # using istarmap, a patch for starmap to allow for tqdm to work
    # get coverage and mutation count one chromosome per thread
    vars_iterable = [(file, chrom_nr, chroms[chrom_nr], clip_filter, clip_filter_3, clip_3, return_read_stats, allow_polyA_clipping) for chrom_nr in range(len(chroms))]
    with Pool(processes=n_cpu) as pool:
        for vec in pool.starmap(coverage_mutations_chrom_PE, vars_iterable):
            vec_chroms.append(vec)

    end = time.time()
    print('time elapsed: %.2f min'%((end-start)/60))
    # unpack the individual vectors
    if not return_read_stats:
        return([i[0] for i in vec_chroms], [i[1] for i in vec_chroms], [i[2] for i in vec_chroms], [i[3] for i in vec_chroms])
    else:
        return([i[0] for i in vec_chroms], [i[1] for i in vec_chroms], [i[2] for i in vec_chroms], [i[3] for i in vec_chroms], [i[4] for i in vec_chroms])

def coverage_mutations_chrom_PE(file, chrom_nr, chrom, clip_filter=9999999, clip_filter_3=True, clip_3=0, return_read_stats=False, allow_polyA_clipping=False):
    """
    Calculate coverage and mutation counts for a single chromosome.

    Args:
        file (str): Path to the BAM file.
        chrom_nr (int): Chromosome number.
        chrom (SeqRecord): Bio.SeqRecord representing the chromosome.
        clip_filter (int): Maximum allowed number of clipped bases.
        clip_filter_3 (bool): Whether to filter based on 3'-end clipping.
        clip_3 (int): Number of bases to remove from the 3'-end.
        return_read_stats (bool): Whether to return read statistics.
        allow_polyA_clipping (bool): Whether to allow poly(A) clipping.

    Returns:
        tuple: A tuple containing a list of vectors (one per chromosome): mutation count plus strand, coverage plus strand, mutation count minus srand, coverage minus strand, stats (optional)
    """

    # open bam file and get chromosome names
    bamfile = pysam.AlignmentFile(file, 'rb')

    chrom_names = bamfile.header.references
    
    # initialize the read_map for the specified chromosome
    vec_mut_p = np.zeros(len(chrom), dtype=int)
    vec_cov_p = np.zeros(len(chrom), dtype=int)
    vec_mut_m = np.zeros(len(chrom), dtype=int)
    vec_cov_m = np.zeros(len(chrom), dtype=int)
    read_stats = []

    dms_dict = {'A':1, 'G':1, 'T':1, 'C':1, '0':0, '?':0, '.':0, '1':1, 'N': 0} 
    
    for read1, read2 in read_pair_generator(bamfile, reference=chrom_names[chrom_nr]):

        if not read1.is_proper_pair or read1.is_secondary or read1.is_supplementary:
            continue

        read_map_1, clipped_3_1, clipped_5_1 = map_mutations(read1, with_clipped=True)
        read_map_2, clipped_3_2, clipped_5_2 = map_mutations(read2, with_clipped=True)
        
        # skip reads with 3'-end clipping >1 (5'-end of read2 is the 3'-end of the molecule)
        if clip_filter_3 and len(clipped_5_2) > 1:
            continue

        # skip reads that have a specified amount of clipped bases elsewhere
        if (len(clipped_3_1) > clip_filter) or (len(clipped_5_1) > clip_filter) or (len(clipped_3_2) > clip_filter) or (len(clipped_5_2) > clip_filter):
            # unless we allow poly(A) clipping in (e.g. mature RNA)
            if allow_polyA_clipping and (is_polyA(clipped_5_2) or is_polyA(clipped_3_1)) and not ((len(clipped_5_1) > clip_filter) or (len(clipped_3_2) > clip_filter)):
                pass
            else:
                continue 

        read_map = join_reads(read_map_1, read_map_2)
        
        n_muts = np.sum(np.isin(np.array(list(read_map.values())), ['A', 'G', 'C', 'T', '1']))
        read_stats.append(n_muts)

        # max. 10% of the read can be mutated
        if n_muts > 0.15 * len(read_map):
            continue

        # PLUS strand
        if not read1.is_reverse:
            # remove specified number of bases from the 3'-end (because Pol protected)
            if clip_3 > 0:
                # find max bp value and remove specified number of previous bases
                end3 = max(read_map.keys())
                [read_map.pop(i, None) for i in np.arange(end3-(clip_3-len(clipped_5_2)), end3+1)]
            for pos in read_map:
                if not read_map[pos] == '.': # don't need to care about clipped bases for now
                    vec_mut_p[pos] += dms_dict[read_map[pos]]
                    vec_cov_p[pos] += 1
        # MINUS strand
        else:
            # remove specified number of bases from the 3'-end (because Pol protected)
            if clip_3 > 0:
                # find max bp value and remove specified number of previous bases
                end3 = min(read_map.keys())
                [read_map.pop(i, None) for i in np.arange(end3, end3+(clip_3-len(clipped_5_2))+1)]
            for pos in read_map:
                if not read_map[pos] == '.': # don't need to care about clipped bases for now
                    vec_mut_m[pos] += dms_dict[read_map[pos]]
                    vec_cov_m[pos] += 1
        
    bamfile.close()

    if not return_read_stats:
        return(vec_mut_p, vec_cov_p, vec_mut_m, vec_cov_m)
    else:
        return(vec_mut_p, vec_cov_p, vec_mut_m, vec_cov_m, read_stats)


def coverage_PE(file, genome_file, n_cpu=4, clip_filter=9999999):
    """
    Parallel wrapper for coverage_chrom_PE.

    Args:
        file (str): Path to the BAM file.
        genome_file (str): Path to the genome file in FASTA format.
        n_cpu (int): Number of CPU cores for parallel processing.
        clip_filter (int): Maximum allowed number of clipped bases.

    Returns:
        tuple: A tuple containing a list of vectors (one per chromosome) for coverage on the plus strand and coverage on the minus strand.
    """

    start = time.time()

    # load the genome    
    chroms = []
    for i in SeqIO.parse(genome_file, 'fasta'):
        chroms.append(i)
    
    vec_chroms = []
    
    # using istarmap, a patch for starmap to allow for tqdm to work
    # get coverage and mutation count one chromosome per thread
    vars_iterable = [(file, chrom_nr, chroms[chrom_nr], clip_filter) for chrom_nr in range(len(chroms))]
    with Pool(processes=n_cpu) as pool:
        for vec in pool.starmap(coverage_chrom_PE, vars_iterable):
            vec_chroms.append(vec)

    end = time.time()
    print('time elapsed: %.2f min'%((end-start)/60))
    # unpack the individual vectors
    return([i[0] for i in vec_chroms], [i[1] for i in vec_chroms])  

def coverage_chrom_PE(file, chrom_nr, chrom, clip_filter=9999999, allow_polyA_clipping=True):
    """
    Calculate coverage for a single chromosome.

    Args:
        file (str): Path to the BAM file.
        chrom_nr (int): Chromosome number.
        chrom (SeqRecord): Bio.SeqRecord representing the chromosome.
        clip_filter (int): Maximum allowed number of clipped bases.
        allow_polyA_clipping (bool): Whether to allow poly(A) clipping.

    Returns:
        tuple: A tuple containing vectors for coverage on the plus strand and coverage on the minus strand.
    """

    # open bam file and get chromosome names
    bamfile = pysam.AlignmentFile(file, 'rb')

    chrom_names = bamfile.header.references
    
    # initialize the read_map for the specified chromosome
    vec_cov_p = np.zeros(len(chrom), dtype=int)
    vec_cov_m = np.zeros(len(chrom), dtype=int)
    
    for read in bamfile.fetch(reference=chrom_names[chrom_nr]):

        read_map, clipped_3, clipped_5 = map_mutations(read, with_clipped=True)
        

        # skip reads that have a specified amount of clipped bases elsewhere
        if (len(clipped_5) > clip_filter) or (len(clipped_3) == 0):
            continue 

        if not is_polyA(clipped_3):
            continue
        
        n_muts = np.sum(np.isin(np.array(list(read_map.values())), ['A', 'G', 'C', 'T', '1']))

        # max. 10% of the read can be mutated
        if n_muts > 0.1 * len(read_map):
            continue

        # PLUS strand
        if not read.is_reverse:
            # remove specified number of bases from the 3'-end (because Pol protected)
            for pos in read_map:
                if not read_map[pos] == '.': # don't need to care about clipped bases for now
                    vec_cov_p[pos] += 1
                    
        # MINUS strand
        else:
            # remove specified number of bases from the 3'-end (because Pol protected)
            for pos in read_map:
                if not read_map[pos] == '.': # don't need to care about clipped bases for now
                    vec_cov_m[pos] += 1
        
    bamfile.close()

    return(vec_cov_p, vec_cov_m)


def join_reads(read1, read2):
    """
    Combine two read maps into a single read map.

    Parameters:
    - read1 (dict): Dictionary representing the first read map.
    - read2 (dict): Dictionary representing the second read map.

    Returns:
    - Dictionary representing the combined read map.
    """
    
    miss_info, ambig_info = '.', '?'
    
    for pos in read2.keys():
        # non-overlapping nts
        if pos not in read1: # add new read2 nts
            read1[pos] = read2[pos]
        # overlapping nts
        else:
            if read1[pos] == read2[pos]:
                continue
            elif read1[pos] not in [miss_info, ambig_info] and read2[pos] in [miss_info, ambig_info]:
                continue
            elif read1[pos] in [miss_info, ambig_info] and read2[pos] not in [miss_info, ambig_info]:
                read1[pos] = read2[pos]
            elif read1[pos] in [miss_info, ambig_info] and read2[pos] in [miss_info, ambig_info]:
                read1[pos] = ambig_info
            else: # disagreement
                read1[pos] = ambig_info
    return(read1)



def map_mutations(read, qscore_cutoff=20, sur_bases=10, with_clipped=False):
    """
    Convert a read's sequence to a vector of 0s & 1s and substituted bases.

    Parameters:
    - read (pysam.AlignedSegment): Read (pysam.AlignedSegment() object).
    - qscore_cutoff (int): Minimum quality score to consider a base.
    - sur_bases (int): Number of surrounding bases to consider for deletion mapping.
    - with_clipped (bool): Whether to return clipped bases as well.

    Returns:
    - Tuple containing vector of 0s & 1s representing mutated positions and substituted bases.
    """

    # 3'-end clipped bases are marked in the read map, 5'-end bases are not
    if read.is_unmapped:
        return('')

    read_map = {}  # Mapping of read to 0s and 1s
    clipped_3 = ''
    clipped_5 = ''

    # if the read is on the - strand we have to get the complementary base
    if read.is_reverse:
        reverse = True
        read_seq = reverse_complement(read.get_forward_sequence()).upper()
    else:
        reverse = False
        read_seq = read.get_forward_sequence().upper()  # Sequence of the read
    ref_seq = read.get_reference_sequence().upper()  # Sequence of the ref genome
    q_scores = read.get_forward_qualities()  # Qual scores of the bases in the read
    ref_pos = fix_ref_pos(read.get_reference_positions(full_length=True))  # Pos in the ref sequence


    i = 0  # Pos in the ref sequence
    j = 0  # Pos in the read sequence
    l = 0  # Pos in the ref position list
    cigar = cigar =re.findall(r'(\d+)([A-Z]{1})', read.cigarstring)
    # print(cigar)
    op_index = 0
    while op_index < len(cigar):  # Each CIGAR operation
        op = cigar[op_index]
        desc, length = op[1], int(op[0])

        if desc == 'M':  # Match or mismatch
            for _ in range(length):  # Each base
                if q_scores[j] >= qscore_cutoff:
                    if not reverse: # register the base that was mutated, depending on strand
                        read_map[ref_pos[l]] = ref_seq[i] if read_seq[j] != ref_seq[i] else '0'
                    else:
                        read_map[ref_pos[l]] = ref_seq[i] if read_seq[j] != ref_seq[i] else '0'
                else:  # < Qscore cutoff
                    read_map[ref_pos[l]] = '?'
                i += 1  # Update ref index
                j += 1  # Update read index
                l += 1  # Update position index

        elif desc == 'D': # Deletion
            if ref_pos[l-1] is None:  # if insertion is followed directly by deletion
                break

            for k in range(length - 1):  # All bases except the 3' end
                read_map[ref_pos[l-1]+k+1] = '?'
                i += 1  # Update ref index
            # 3' end of deletion
            ambig = map_deletion(ref_seq, i, length, sur_bases)
            read_map[ref_pos[l-1]+length] = '?' if ambig else '1'
            i += 1  # Update ref index

        elif desc == 'I':  # Insertion
            j += length  # Update read index
            l += length  

        elif desc == 'S':  # Soft clipping
            if (not reverse and op_index == len(cigar) - 1) or (reverse and op_index == 0):  # Soft clipped at the 3'-end
                for _ in range(length):
                    read_map[ref_pos[l]] = '.'
                    l += 1  # Update position index (soft-clipped bases are part of the position list)
                clipped_3 = read_seq[j:j+length] # clipped 3'-end bases for studying poly(A) tail properties
            else: # Soft clipped at 5'-end
                for _ in range(length):
                    read_map[ref_pos[l]] = '.'
                    l += 1  # Update position index (soft-clipped bases are part of the position list)
                clipped_5 = read_seq[j:j+length] # clipped 5'-end bases for studying RT-stop tail properties
            j += length  # Update read index

        elif desc == 'N': # Intron, already accounted for by pysam in the ref seq
            pass

        elif desc == 'H': # hard clip, already accounted for by minimap2
            pass

        else:
            print('Unknown CIGAR op encountered: %s'%(desc))
            break

        op_index += 1

    # return the vector and clipped bases if desired
    if with_clipped:
        if reverse:
            return(read_map, reverse_complement(clipped_3), reverse_complement(clipped_5))
        else:
            return(read_map, clipped_3, clipped_5)
    else:
        return(read_map)

def map_deletion(ref_seq, i, length, num_surBases):
    """
    Determines whether a deletion is ambiguous or not by looking at the sequence surrounding the deletion.
    
    Parameters:
        ref_seq (str): Reference sequence
        i (int): 3' index of deletion in reference sequence
        length (int): Length of deletion
        num_surBases (int): Number of surrounding bases to consider
    
    Returns:
        bool: Whether deletion is ambiguous or not
    """

    orig_del_start = i - length + 1
    orig_sur_start = orig_del_start - num_surBases
    orig_sur_end = i + num_surBases
    orig_sur_seq = ref_seq[orig_sur_start - 1: orig_del_start - 1] + ref_seq[i:orig_sur_end]
    for new_del_end in range(i - length, i + length + 1):  # Alt del end points
        if new_del_end == i:  # Orig end point
            continue
        new_del_start = new_del_end - length + 1
        sur_seq = ref_seq[orig_sur_start - 1: new_del_start - 1] + ref_seq[new_del_end:orig_sur_end]
        if sur_seq == orig_sur_seq:
            return True
    return False


def fix_ref_pos(ref_pos):
    """
    Fills in None values at the soft-clipped positions that pysam produces.
    
    Parameters:
        ref_pos (list): List of soft-clipped positions.
    
    Returns:
        list: Updated list with filled soft-clipped positions.
    """

    if None not in ref_pos:
        return(ref_pos)

    i = 0
    clip_counter = 0
    while ref_pos[i] is None:
        i += 1
        clip_counter += 1
    
    if i > 0:
        for j in range(clip_counter):
            ref_pos[j] = ref_pos[i] - clip_counter + j
   
    if i == len(ref_pos)-1:
        return(ref_pos)
    
    ref_pos = ref_pos[::-1]
    i = 0
    clip_counter = 0
    while ref_pos[i] is None:
        i += 1
        clip_counter += 1

    if i > 0:
        for j in range(clip_counter):
            ref_pos[j] = ref_pos[i] + clip_counter - j
    
    return(ref_pos[::-1])

def auroc(x, y):  
    """
    Calculate AUROC (Area Under the Receiver Operating Characteristic curve) using the U statistic.
    
    Parameters:
        x (array-like): Data for the first group.
        y (array-like): Data for the second group.
    
    Returns:
        float: AUROC value.
    """

    mwu = mannwhitneyu(x, y)
    u = mwu[0]
    n1, n2 = len(x), len(y)
    area = u/(n1*n2)
    return(area)

def gini(arr):
    """
    Calculate the Gini coefficient of a numeric array.
    
    Parameters:
        arr (array-like): Numeric array.
    
    Returns:
        float: Gini coefficient.
    """

    if len(arr) == 0:
        return(np.nan)

    accum, giniB = 0, 0
    for i in sorted(arr):
        accum += i
        giniB += accum - i / 2.0
    fair_area = accum * len(arr) / 2.0
    return((fair_area - giniB) / fair_area)


def is_polyA(seq, min_A_frac=0.8):
    """
    Check if a given sequence is a polyadenine (polyA) sequence.

    Parameters:
        seq (str): Input sequence.
        min_A_frac (float): Minimum fraction of 'A' bases for considering it a polyA sequence.

    Returns:
        bool: True if it's a polyA sequence, False otherwise.
    """

    c = Counter(seq)
    if c and c['A'] > 0 and c['A']/sum(c.values()) >= min_A_frac:
        return(True)
    else:
        return(False)


def rsample(seq, rea, msk, rnastructure_path, md=9999, ns=10000, t=303.15, max_rea=1):
    """
    Perform RNA structure prediction using sampling and return shannon entropy, free energy, and predicted structure.

    Parameters:
        seq (str): RNA sequence.
        rea (array-like): Reactivities.
        msk (array-like): Mask for reactivities.
        rnastructure_path (str): Path to RNAstructure executable.
        md (int): Maximum pairing distance for RNAstructure.
        ns (int): Number of samples for RNAstructure.
        t (float): Temperature in K.
        max_rea (int): Maximum reactivity value.

    Returns:
        tuple: Tuple containing shannon entropy, free energy, and predicted structure.
    """

    fd_shape, path_shape = tempfile.mkstemp()
    fd_seq,   path_seq   = tempfile.mkstemp()
    fd_out,   path_out   = tempfile.mkstemp()
    fd_prob,  path_prob  = tempfile.mkstemp()
    fd_ct,    path_ct    = tempfile.mkstemp()
    fd_dbr,   path_dbr   = tempfile.mkstemp()

    try: 
        
        # write the reactivities for RNAstructure to temp file
        with os.fdopen(fd_shape, 'w') as f:
            for i in range(len(rea)):
                if msk[i]:
                    f.write(str(i+1) + '\t' + str(np.min([rea[i], max_rea])) + '\n')
                else:
                    f.write(str(i+1) + '\t' + str(-999) + '\n')

        # write the sequence for RNAstructure to temp file
        with os.fdopen(fd_seq, 'w') as f:
            f.write('>generic' + '\n')
            f.write(seq)

        # calculate partition function
        out = subprocess.run([rnastructure_path + '/Rsample-smp', path_seq, path_shape, path_out, '-md', str(md), '-ns', str(ns), '-t', str(t), '--DMS'], capture_output=False, shell=False, stdout=open(os.devnull, 'wb'))

        # calculate MEA
        out = subprocess.run([rnastructure_path + '/EnsembleEnergy', path_out], capture_output=True, shell=False)
        E = float(out.stdout.decode().split('\n')[3].split(' ')[-2])
        
        # calculate pairing probabilities
        out = subprocess.run([rnastructure_path + '/ProbabilityPlot', path_out, path_prob, '-t'], capture_output=False, shell=False, stdout=open(os.devnull, 'wb'))
        probs = pd.read_csv(path_prob, delimiter='\t', skiprows=1)

        # calculate most likely pairings
        out = subprocess.run([rnastructure_path + '/MaxExpect-smp', path_out, path_ct], capture_output=False, shell=False, stdout=open(os.devnull, 'wb'))
        out = subprocess.run([rnastructure_path + '/ct2dot', path_ct, '1', path_dbr], capture_output=False, shell=False, stdout=open(os.devnull, 'wb'))
        with open(path_dbr, 'r') as f:
            f.readline()
            f.readline()
            dbr = f.readline().strip()

        # calculate shannon entropy for each nucleotide
        shan = np.zeros(len(seq))
        done = []
        for idx, pair in probs.iterrows():
            shan[int(pair['i'])-1] += 10**(-pair['-log10(Probability)']) * -pair['-log10(Probability)']
            shan[int(pair['j'])-1] += 10**(-pair['-log10(Probability)']) * -pair['-log10(Probability)']
            done.append(pair['i']-1)
            done.append(pair['j']-1)
        
        shan = -np.array(shan)
        shan[~np.isin(np.arange(0, len(seq)), done)] = np.nan 
        # print(shan)
        # print(E)

    finally:
        os.remove(path_shape)
        os.remove(path_seq)
        os.remove(path_out)
        os.remove(path_prob)
        os.remove(path_ct)
        os.remove(path_dbr)

    return(shan, E, dbr)


def fold_dms(seq, rea, msk, rnastructure_path, md=9999, m=3, t=303.15, l=30, max_rea=3):
    """
    Perform RNA structure folding with DMS reactivity data.

    Parameters:
        seq (str): RNA sequence.
        rea (array-like): Reactivities.
        msk (array-like): Mask for reactivities.
        rnastructure_path (str): Path to RNAstructure executable.
        md (int): Maximum pairing distance for RNAstructure.
        m (int): Number of suboptimal structures to generate.
        t (float): Temperature in K.
        l (int): Maximum base pair span for DMS reactivity.
        max_rea (int): Maximum reactivity value.

    Returns:
        list: List of suboptimal structures in dot-bracket notation.
    """

    fd_shape, path_shape = tempfile.mkstemp()
    fd_seq,   path_seq   = tempfile.mkstemp()
    fd_ct,    path_ct    = tempfile.mkstemp()
    fd_dbr,   path_dbr   = tempfile.mkstemp()

    try: 
        
        # write the reactivities for RNAstructure to temp file
        with os.fdopen(fd_shape, 'w') as f:
            for i in range(len(rea)):
                if msk[i]:
                    f.write(str(i+1) + '\t' + str(np.min([rea[i], max_rea])) + '\n')
                else:
                    f.write(str(i+1) + '\t' + str(-999) + '\n')

        # write the sequence for RNAstructure to temp file
        with os.fdopen(fd_seq, 'w') as f:
            f.write('>generic' + '\n')
            f.write(seq)

        # call fold
        out = subprocess.run([rnastructure_path + '/Fold-smp', path_seq, path_ct, '-dms', path_shape, '-md', str(md), '-m', str(m), '-t', str(t), '-l', str(l)], capture_output=False, shell=False, stdout=open(os.devnull, 'wb'))
        
        # convert to dbr
        dbrs = []
        for i in range(m):
            out = subprocess.run([rnastructure_path + '/ct2dot', path_ct, str(i+1), path_dbr], capture_output=False, shell=False, stdout=open(os.devnull, 'wb'))
            with open(path_dbr, 'r') as f:
                f.readline()
                f.readline()
                dbrs.append(f.readline().strip())


    finally:
        os.remove(path_shape)
        os.remove(path_seq)
        os.remove(path_ct)
        os.remove(path_dbr)

    return(dbrs)


def fold(seq, t=303.15, maxbp=200, rnastructure_path='/Users/leo/Builds/RNAstructure/exe'):
    """
    Fold an RNA sequence and return the dot-bracket notation.

    Parameters:
        seq (str): RNA sequence.
        t (float): Temperature in K.
        maxbp (int): Maximum base pair span.
        rnastructure_path (str): Path to RNAstructure executable.

    Returns:
        tuple: Tuple containing dot-bracket notation and free energy.
    """

    # call folding function
    out = subprocess.run([rnastructure_path + '/Fold', '-', '-', '-t', str(t), '--bracket', '--mfe', '-md', str(maxbp)], capture_output=True, shell=False, input=f">shart\n{seq.upper()}".encode())
    
    # read output
    out = out.stdout.decode().split('\n')

    if len(out[0].split(' ')) == 1:
        E = 999
    else:
        E = float(out[0].split(' ')[2])

    dbr = out[2]

    return(dbr, E)

def fold_parallel(seqs, t=303.15, maxbp=200, rnastructure_path='/Users/leo/Builds/RNAstructure/exe', n_cpu=6):
    """
    Fold multiple RNA sequences in parallel.

    Parameters:
        seqs (list): List of RNA sequences.
        t (float): Temperature in K.
        maxbp (int): Maximum base pair span.
        rnastructure_path (str): Path to RNAstructure executable.
        n_cpu (int): Number of CPU cores to use.

    Returns:
        tuple: Tuple containing lists of dot-bracket notations and free energies.
    """

    start = time.time()
    
    dbr_all = []
    eng_all = []

    vars_iterable = [(seq, t, maxbp, rnastructure_path) for seq in seqs]
    with Pool(processes=n_cpu) as pool:
        for dbr, e in pool.starmap(fold, vars_iterable):
            dbr_all.append(dbr)
            eng_all.append(e)

    end = time.time()
    print('time elapsed: %.2f min'%((end-start)/60))
    return(dbr_all, eng_all)


def dynamic_auroc(rea, dbr, mask_seq, n_min, offset):
    """
    Calculate AUROC in dynamic windows based on RNA reactivity and structure.

    Parameters:
        rea (array-like): Reactivities.
        dbr (str): Dot-bracket notation.
        mask_seq (array-like): Mask for sequence positions.
        n_min (int): Minimum number of nucleotides in each window.
        offset (int): Offset for window scanning.

    Returns:
        tuple: Tuple containing lists of AUROC values, window indices, and overall AUROC.
    """
    
    # finds all possible windows of adjusted input size to calculate the auroc in
    np.seterr(all="ignore")
    mask_un = np.array([True if i == '.' else False for i in dbr]) & mask_seq
    mask_pa = np.array([True if i in [')', '('] else False for i in dbr]) & mask_seq
    mask_al = mask_pa | mask_un
    indices = np.arange(0, len(dbr))[mask_al]
    mask_fi = mask_pa[mask_al]
    
    # only consider A/C and bases with structure information
    avg = rea[mask_al]
    
    # find optimal window composition by comparing total number of paired/unpaired bases
    if np.sum(mask_un) > np.sum(mask_pa):
        n_pa_target = n_min
        n_un_target = np.round(n_min * (np.sum(mask_un) / np.sum(mask_pa)))
    elif np.sum(mask_un) < np.sum(mask_pa):
        n_un_target = n_min
        n_pa_target = np.round(n_min * (np.sum(mask_pa) / np.sum(mask_un)))
    else:
        n_pa_target = n_min
        n_un_target = n_min

    i = offset
    auroc_all = []
    index_all = []
    done      = False
    while not done:
        # find the next window that satisfies precomputed number of nts
        rea_un = []
        rea_pa = []
        idx_al = []
        while len(rea_pa) < n_pa_target or len(rea_un) < n_un_target:
            if i >= len(mask_fi)-1: # stop if it happens to perfectly hit the end (doesn't happen often)
                break
                
            if mask_fi[i]:
                rea_pa.append(avg[i])
            else:
                rea_un.append(avg[i])
            idx_al.append(indices[i])
            i += 1

        # calculate the auroc and check if conditions are met
        if len(rea_pa) >= n_pa_target and len(rea_un) >= n_un_target: # only in case the break statement is executed
            auroc_all.append(auroc(rea_un, rea_pa))
            index_all.append(idx_al)

        if len(avg[i:]) < n_pa_target + n_un_target:
            done = True
    
    # print(f"Target size for windows is {n_pa_target:.0f} paired and {n_un_target:.0f} unpaired nts.")
    # print(f"Found {len(auroc_all)} windows with average size {np.mean([len(i) for i in index_all]):.2f}")
    auroc_overall = auroc(avg[~mask_fi], avg[mask_fi])
    
    return(auroc_all, index_all, auroc_overall)

def auroc_all_windows(rea, dbr, mask_seq, n_min, n_cpu=5):
    """
    Calculate AUROC in dynamic windows for all possible offsets in parallel.

    Parameters:
        rea (array-like): Reactivities.
        dbr (str): Dot-bracket notation.
        mask_seq (array-like): Mask for sequence positions.
        n_min (int): Minimum number of nucleotides in each window.
        n_cpu (int): Number of CPU cores to use.

    Returns:
        tuple: Tuple containing lists of window indices, AUROC values, and overall AUROC.
    """

    start = time.time()
    
    auroc_all = []
    index_all = []
    # calculate auroc in windows for this particular offset, one per thread
    vars_iterable = [(rea, dbr, mask_seq, n_min, i) for i in np.arange(0, dbr.count('.')+dbr.count('(')+dbr.count(')'), 1)]
    with Pool(processes=n_cpu) as pool:
        for auroc_off, index_off, auroc_overall in pool.starmap(dynamic_auroc, vars_iterable):
            # check if this window has been found before
            # this is brute force and therefore slow, whatever. what you gonna do. it's parallel so f*** off
            for aur, ind in zip(auroc_off, index_off):
                if not ind in index_all:
                    index_all.append(ind)
                    auroc_all.append(aur)

    end = time.time()
    print('time elapsed: %.2f min'%((end-start)/60))
    return(index_all, auroc_all, auroc_overall)


def sliding_corrcoeff(rea1, rea2, w):
    """
    Calculate sliding window correlation coefficients between two reactivity arrays.

    Parameters:
        rea1 (array-like): Reactivity array 1.
        rea2 (array-like): Reactivity array 2.
        w (int): Window size.

    Returns:
        list: List of correlation coefficients.
    """

    # w must be even number so nt can be in the middle
    if w%2 != 0:
        w += 1
        print(f"Adjusted window size to {w}.")
        
    assert len(rea1) == len(rea2), "Input arrays must have the same length."
         
    r_all = []
    
    for i in np.arange(0, len(rea1), 1):
        if i < w/2:
            r_all.append(np.nan)
        elif i >= len(rea1)-1-w/2:
            r_all.append(np.nan)
        else:
            r_all.append(pearsonr(rea1[i-int(w/2):i+int(w/2)], rea2[i-int(w/2):i+int(w/2)]).statistic)
        
    return(r_all)


def normalize_cov(mut, cov, return_params=False):
    """
    Normalize mutation and coverage data using the method of moments.

    Parameters:
        mut (array-like): Mutation data.
        cov (array-like): Coverage data.
        return_params (bool): Whether to return alpha_p and beta_p.

    Returns:
        array-like or tuple: Normalized mutation data or tuple of alpha_p and beta_p.
    """

    # Estiread alpha_p and beta_p using method of moments
    rat = (mut + 1) / (cov + 1)
    E = np.mean(rat)
    V = np.var(rat)
    alpha_p = E * (((E * (1 - E)) / V) - 1)
    beta_p = (alpha_p * (1 - E)) / E

    # Estiread adjusted mutation rate
    rat_adj = (mut + alpha_p) / (cov + beta_p)

    if not return_params:
        return(rat_adj)
    else:
        return(alpha_p, beta_p)


def weighted_moving_average(t, y, window_size):
    """
    Calculate a weighted moving average.

    Parameters:
        t (array-like): "Time" array.
        y (array-like): Value array.
        window_size (int): Size of the moving window.

    Returns:
        array-like: Weighted moving averages.
    """

    if len(t) != len(y):
        raise ValueError("Time and value arrays must have the same length.")

    n = len(y)
    averages = [0] * n

    for i in range(n):
        total_weight = 0
        weighted_sum = 0

        # Iterate from the current index `i` to a maximum of `max(-1, i - window_size)`
        for j in range(i, max(-1, i - window_size), -1):
            if y[j] is not None:
                weight = t[i] - t[j]
                weighted_sum += y[j] * weight
                total_weight += weight

        if total_weight != 0:
            averages[i] = weighted_sum / total_weight

    return(np.array(averages))


def moving_average(y, window_size):
    """
    Calculate a simple moving average.

    Parameters:
        y (array-like): Value array.
        window_size (int): Size of the moving window.

    Returns:
        array-like: Simple moving averages.
    """

    avgs = []

    assert window_size % 2 == 0

    for i in np.arange(0.5*window_size, len(y)-0.5*window_size+1):
        window_data = y[int(i)-int(0.5*window_size):int(i)+int(0.5*window_size)]
        if np.sum(~np.isnan(window_data)) >= 3:
            avgs.append(np.nanmean(window_data))
        else:
            avgs.append(np.nan)

    head = np.array([avgs[0]]*int(0.5*window_size))
    tail = np.array([avgs[-1]]*(int(0.5*window_size)-1))

    return(np.hstack([head, avgs, tail]))


def get_coords(seq, dbr, vienna_location='/Users/leo/Builds/ViennaRNA/bin'):
    """
    Get coordinates from RNA secondary structure prediction using ViennaRNA with the RNA Puzzler algorithm.

    Parameters:
        seq (str): RNA sequence.
        dbr (str): Dot-bracket notation.
        vienna_location (str): Path to ViennaRNA executable.

    Returns:
        array-like: Puzzler coordinates.
    """

    assert len(seq) == len(dbr), 'Input variables must have same length.'
    
    fd_seq, path_seq = tempfile.mkstemp()
    
    try:
        # Export the prediction for ViennaRNA
        with os.fdopen(fd_seq, 'w') as f:
            f.write('>temp_rna file\n')
            f.write(seq + '\n')
            f.write(dbr + '\n')

        # RNAplot from ViennaRNA is called here to get the puzzler coordinates
        exit = os.system(f"cat {path_seq} | {vienna_location}/RNAplot -t 4")
        if exit != 0:
            raise(Exception('Bracket notation not coherent.'))

        # parse the postscript output file and collect the coordinates
        puzzler = []
        with open('temp_rna_ss.ps', 'r') as f:
            found = False
            over = False
            line = f.readline()
            while not found:
                found = line == '/coor [\n'
                line = f.readline()
            while not over:
                puzzler.append(line.strip().strip('[]').split(' '))
                line = f.readline()
                over = line == '] def\n'
                
    finally:
        os.remove(path_seq)
        os.remove('temp_rna_ss.ps')

    return(np.array(puzzler, dtype=float))

def plot_structure(coords, seq, rea, mask, axs=None, cmap_loc='cmap.txt', circle_size=50, text_size=7, figsize=(10,8), vmin=0, vmax=1):
    """
    Plot RNA secondary structure with reactivity data.

    Parameters:
        coords (array-like): Puzzler coordinates.
        seq (str): RNA sequence.
        rea (array-like): Reactivities.
        mask (array-like): Mask for sequence positions.
        axs (matplotlib.axes.Axes): Matplotlib axes.
        cmap_loc (str): Path to colormap file.
        circle_size (int): Size of circles in the plot.
        text_size (int): Size of text labels in the plot.
        figsize (tuple): Figure size.
        vmin (float): Minimum value for colormap.
        vmax (float): Maximum value for colormap.
    """
    
    if not axs:
        fig, axs = plt.subplots(figsize=figsize)
        fig.patch.set_facecolor('white')

    # get colors   
    cmap = np.array(pd.read_csv(cmap_loc, header=None))/255
    cmap = mplc.ListedColormap(cmap)
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    m = cm.ScalarMappable(norm=norm, cmap=cmap)
    clrs = m.to_rgba(rea)
    clrs[~mask] = [230/255, 230/255, 230/255, 1]

    axs.plot(coords[:,0], coords[:,1], '-k', zorder=0)
    axs.scatter(coords[:,0], coords[:,1], circle_size, color=clrs, edgecolors='none')
    for i, nt in enumerate(seq.replace('T', 'U')):
        axs.text(coords[i,0], coords[i,1], nt, ha='center', va='center', color='k', size=text_size)
    axs.axis('equal')
    axs.axis('off')






