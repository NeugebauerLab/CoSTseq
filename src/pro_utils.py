from src.bam_utils import read_pair_generator, read_ends, get_clipped
from src.dms_utils import map_mutations, join_reads
import pysam
import numpy as np
from collections import Counter
from multiprocessing import Pool
from Bio import SeqIO, Seq
import pandas as pd
import time
from tqdm import tqdm
from scipy.stats import binned_statistic

def DMS_read2(file, length, genome_file, chrom_name=None, clip_filter=5):
    """
    Calculate the mutation and coverage vectors only from read2 alignments.

    Parameters:
    - file (str): Path to the BAM file containing alignments.
    - length (int): Length of the genomic region.
    - genome_file (str): Path to the FASTA file containing the reference genome.
    - chrom_name (str): Name of the chromosome to analyze.
    - clip_filter (int): Maximum allowed number of clipped bases.

    Returns:
    - tuple: Two 1D NumPy arrays representing mutation and coverage vectors for the specified.
    """

    # load chromosome sequence
    for c in SeqIO.parse(genome_file, 'fasta'):
        if c.id == chrom_name:
            seq = str(c.seq)

    bamfile = pysam.AlignmentFile(file, 'rb')

    vec_cov = np.zeros(length, dtype=int)
    vec_mut = np.zeros(length, dtype=int)

    for read in bamfile.fetch(reference=chrom_name):
        if read.is_read2 and read.is_proper_pair and not read.is_secondary and not read.is_supplementary:

            read_map, clipped_3, clipped_5 = map_mutations(read, with_clipped=True)
            
            # skip junction reads
            if not np.all(np.diff(list(read_map.keys()))==1):
                continue

            # skip reads with 3'-end clipping >1 (5'-end of read2 is the 3'-end of the molecule)
            if len(clipped_5) > 1:
                continue

            # skip reads that have a specified amount of clipped bases at 5'
            if len(clipped_3) > clip_filter:
                continue

            # max. 10% of the read can be mutated
            if np.sum(np.isin(np.array(list(read_map.values())), ['A', 'G', 'C', 'T', '1'])) > 0.1 * len(read_map):
                continue

            if read.is_reverse:
                for i, pos in enumerate(reversed(read_map)):
                    pos_norm = list(read_map.keys())[-1]-pos
                    if pos_norm >= length:
                        break
                    
                    if read_map[pos] == '0' and seq[pos] in ['A', 'C']:
                        vec_cov[pos_norm] += 1
                            
                    elif (read_map[pos] in ['A', 'C']) or (read_map[pos] == '1' and seq[pos] in ['A', 'C']):
                        vec_mut[pos_norm] += 1
                        vec_cov[pos_norm] += 1

                    elif i == 0 and read_map[pos] == '.' and seq[pos] in ['A', 'C']: # the one allowed clipped base should be registered as a mutation
                        vec_mut[pos_norm] += 1
                        vec_cov[pos_norm] += 1                
                    
            else:
                for i, pos in enumerate(read_map):
                    pos_norm = pos-list(read_map.keys())[0]
                    if pos_norm >= length:
                        break
                    
                    if read_map[pos] == '0' and seq[pos] in ['G', 'T']:
                        vec_cov[pos_norm] += 1
                            
                    elif (read_map[pos] in ['G', 'T']) or (read_map[pos] == '1' and seq[pos] in ['G', 'T']):
                        vec_mut[pos_norm] += 1
                        vec_cov[pos_norm] += 1

                    elif i == 0 and read_map[pos] == '.' and seq[pos] in ['G', 'T']: # the one allowed clipped base should be registered as a mutation
                        vec_mut[pos_norm] += 1
                        vec_cov[pos_norm] += 1 
    
    bamfile.close()    
    return(vec_mut, vec_cov)

def DMS_read2_parallel(file, length, genome_file, n_cpu=4, clip_filter=5):
    """
    Parallel version of DMS_read2 function.

    Parameters:
    - file (str): Path to the BAM file containing read2 alignments.
    - length (int): Length of the genomic region.
    - genome_file (str): Path to the FASTA file containing the reference genome.
    - n_cpu (int): Number of CPU cores to use.
    - clip_filter (int): Maximum allowed number of clipped bases at the 5' end.

    Returns:
    - tuple: Two 1D NumPy arrays representing mutation and coverage vectors.
    """
    
    start = time.time()

    bamfile = pysam.AlignmentFile(file, 'rb')
    chrom_names = bamfile.header.references
    bamfile.close()

    vec_cov_all = np.zeros(length, dtype=int)
    vec_mut_all = np.zeros(length, dtype=int)

    vars_iterable = [(file, length, genome_file, chrom_name, clip_filter) for chrom_name in chrom_names]
    with Pool(processes=n_cpu) as pool:
        for vec_mut, vec_cov in pool.starmap(DMS_read2, vars_iterable):
            vec_cov_all += vec_cov
            vec_mut_all += vec_mut

    end = time.time()
    print('time elapsed: %.2f min'%((end-start)/60))

    return(vec_mut_all, vec_cov_all)


def DMS_read2_region(filename, length, ref, start, end, strand, genome_file, clip_filter=5):
    """
    Calculate the mutation and coverage vectors for a specific genomic region.

    Parameters:
    - filename (str): Path to the BAM file containing alignments.
    - length (int): Length of the genomic region.
    - ref (str): Reference chromosome name.
    - start (int): Start position of the region.
    - end (int): End position of the region.
    - strand (str): Strand of interest ('+' or '-').
    - genome_file (str): Path to the FASTA file containing the reference genome.
    - clip_filter (int): Maximum allowed number of clipped bases.

    Returns:
    - tuple: Two 1D NumPy arrays representing mutation and coverage vectors for the specified genomic region.
    """

    # one-based coordinates
    # load chromosome sequence
    for c in SeqIO.parse(genome_file, 'fasta'):
        if c.id == ref:
            seq = str(c.seq)

    bamfile = pysam.AlignmentFile(filename, 'rb')

    vec_cov = np.zeros(length, dtype=int)
    vec_mut = np.zeros(length, dtype=int)

    try:        
        for read in bamfile.fetch(reference=ref, start=start+1, stop=end+1):
            if read.is_read2 and read.is_proper_pair and not read.is_secondary and not read.is_supplementary:
            
                # strand check
                if (not read.is_reverse and strand == '+') or (read.is_reverse and strand == '-'):
                    continue

                # determine RNA 3'-end position, which is always the 5'-end of read2
                _, mapped_5, mapped_3 = read_ends(read)
                threeprime = mapped_5-start

                # 3'-end outside of requested region 
                if threeprime < 0 or end-mapped_5 <= 0: 
                    continue
                    
                # obtain mutational vector, clipped bases
                read_map, clipped_3, clipped_5 = map_mutations(read, with_clipped=True)

                # skip junction reads
                if not np.all(np.diff(list(read_map.keys()))==1):
                    continue

                # skip reads with 3'-end clipping >1 (5'-end of read2 is the 3'-end of the molecule)
                if len(clipped_5) > 1:
                    continue

                # skip reads that have a specified amount of clipped bases at 5'
                if len(clipped_3) > clip_filter:
                    continue

                # max. 10% of the read can be mutated
                if np.sum(np.isin(np.array(list(read_map.values())), ['A', 'G', 'C', 'T', '1'])) > 0.1 * len(read_map):
                    continue

                if read.is_reverse:
                    for i, pos in enumerate(reversed(read_map)):
                        pos_norm = list(read_map.keys())[-1]-pos
                        if pos_norm >= length:
                            break
                        
                        if read_map[pos] == '0' and seq[pos] in ['A', 'C']:
                            vec_cov[pos_norm] += 1
                                
                        elif (read_map[pos] in ['A', 'C']) or (read_map[pos] == '1' and seq[pos] in ['A', 'C']):
                            vec_mut[pos_norm] += 1
                            vec_cov[pos_norm] += 1

                        elif i == 0 and read_map[pos] == '.' and seq[pos] in ['A', 'C']: # the one allowed clipped base should be registered as a mutation
                            vec_mut[pos_norm] += 1
                            vec_cov[pos_norm] += 1                
                        
                else:
                    for i, pos in enumerate(read_map):
                        pos_norm = pos-list(read_map.keys())[0]
                        if pos_norm >= length:
                            break
                        
                        if read_map[pos] == '0' and seq[pos] in ['G', 'T']:
                            vec_cov[pos_norm] += 1
                                
                        elif (read_map[pos] in ['G', 'T']) or (read_map[pos] == '1' and seq[pos] in ['G', 'T']):
                            vec_mut[pos_norm] += 1
                            vec_cov[pos_norm] += 1

                        elif i == 0 and read_map[pos] == '.' and seq[pos] in ['G', 'T']: # the one allowed clipped base should be registered as a mutation
                            vec_mut[pos_norm] += 1
                            vec_cov[pos_norm] += 1 

    finally:
        bamfile.close()

    return(vec_mut, vec_cov)

def DMS_read2_region_parallel(filename, length, genome_file, regions_csv, n_cpu=4, clip_filter=5):
    """
    Parallel version of DMS_read2_region function.

    Parameters:
    - filename (str): Path to the BAM file containing alignments.
    - length (int): Length of the genomic region.
    - genome_file (str): Path to the FASTA file containing the reference genome.
    - regions_csv (str): Path to the CSV file containing region information.
    - n_cpu (int): Number of CPU cores to use.
    - clip_filter (int): Maximum allowed number of clipped bases at the 5' end.

    Returns:
    - tuple: Two 1D NumPy arrays representing mutation and coverage vectors for the specified genomic region.
    """
    
    start = time.time()

    vec_cov_all = np.zeros(length, dtype=int)
    vec_mut_all = np.zeros(length, dtype=int)

    regions = pd.read_csv(regions_csv, delimiter=',', header=None)

    vars_iterable = [(filename, length, ref, start, end, strand, genome_file, clip_filter) for _, (ref, start, end, strand) in regions.iterrows()]
    with Pool(processes=n_cpu) as pool:
        for vec_mut, vec_cov in pool.starmap(DMS_read2_region, vars_iterable):
            vec_cov_all += vec_cov
            vec_mut_all += vec_mut

    end = time.time()
    print('time elapsed: %.2f min'%((end-start)/60))

    return(vec_mut_all, vec_cov_all)


def map_read_ends(file, genome_file, with_lengths=False, clip_filter=5):
    """
    Map the read ends and calculate statistics.

    Parameters:
    - file (str): Path to the BAM file containing alignments.
    - genome_file (str): Path to the FASTA file containing the reference genome.
    - with_lengths (bool): Whether to include insert lengths in the output.
    - clip_filter (int): Maximum allowed number of clipped bases.

    Returns:
    - tuple: Depending on 'with_lengths', returns various NumPy arrays containing read end statistics.
    """
    
    # load the genome    
    chroms = []
    for i in SeqIO.parse(genome_file, 'fasta'):
        chroms.append(i)
    
    bamfile = pysam.AlignmentFile(file, 'rb')
    chrom_names = bamfile.header.references

    # initialize the bitvector for the specified chromosome
    vec_polym_plus = [np.zeros(len(chrom), dtype=int) for chrom in chroms]
    vec_rstop_plus = [np.zeros(len(chrom), dtype=int) for chrom in chroms]
    vec_polym_minus = [np.zeros(len(chrom), dtype=int) for chrom in chroms]
    vec_rstop_minus = [np.zeros(len(chrom), dtype=int) for chrom in chroms]
    insert_lengths = []

    start = time.time()
     
    for read1, read2 in read_pair_generator(bamfile):

        clipped_3_1, clipped_5_1 = get_clipped(read1)
        clipped_3_2, clipped_5_2 = get_clipped(read2)

        # skip reads with 3'-end clipping >1 (5'-end of read2 is the 3'-end of the molecule)
        if len(clipped_5_2) > 1:
            continue

        # skip reads that have a specified amount of clipped bases elsewhere
        if (len(clipped_3_1) > clip_filter) or (len(clipped_5_1) > clip_filter) or (len(clipped_3_2) > clip_filter):
            continue 

        ref, fivep_r1, _ = read_ends(read1)
        _, fivep_r2, _ = read_ends(read2)

        # if there is on base clipped at 3'-end, shift the position accordingly
        if len(clipped_5_2) == 1 and not read1.is_reverse:
            fivep_r2 += 1
        elif len(clipped_5_2) == 1 and read1.is_reverse:
            fivep_r2 -= 1

        insert_lengths.append(abs(fivep_r1-fivep_r2))

        assert read1.is_read1

        # PLUS strand
        if not read1.is_reverse:
            vec_rstop_plus[chrom_names.index(ref)][fivep_r1] += 1
            vec_polym_plus[chrom_names.index(ref)][fivep_r2] += 1
                
        # MINUS strand
        else:
            vec_rstop_minus[chrom_names.index(ref)][fivep_r1] += 1
            vec_polym_minus[chrom_names.index(ref)][fivep_r2] += 1

    bamfile.close()
    end = time.time()
    print('time elapsed: %.2f min'%((end-start)/60))

    if with_lengths:
        return(vec_polym_plus, vec_rstop_plus, vec_polym_minus, vec_rstop_minus, insert_lengths)
    else:
        return(vec_polym_plus, vec_rstop_plus, vec_polym_minus, vec_rstop_minus)



def three_base_identity(polym_p, polym_m, seq, min_cov=30, surr=5):
    """
    Calculate 3'-end base identity based on polymerase occupancy data.

    Parameters:
    - polym_p (numpy.ndarray): Polymerase occupancy data for the plus strand.
    - polym_m (numpy.ndarray): Polymerase occupancy data for the minus strand.
    - seq (str): RNA sequence.
    - min_cov (int): Minimum coverage threshold.
    - surr (int): Surrounded region for calculation.

    Returns:
    - pandas.DataFrame: DataFrame with normalized frequencies for each base.
    """

    def _shift(arr, num, fill_value=np.nan):
        result = np.empty_like(arr)
        if num > 0:
            result[:num] = fill_value
            result[num:] = arr[:-num]
        elif num < 0:
            result[num:] = fill_value
            result[:num] = arr[-num:]
        else:
            result[:] = arr
        return(result)

    A_freq = [0 for i in range(2*surr+1)]
    G_freq = [0 for i in range(2*surr+1)]
    U_freq = [0 for i in range(2*surr+1)]
    C_freq = [0 for i in range(2*surr+1)]
    for i in range(len(seq)):
        
        A = np.array(list(seq[i])) == 'A'
        G = np.array(list(seq[i])) == 'G'
        U = np.array(list(seq[i])) == 'T'
        C = np.array(list(seq[i])) == 'C'
        Sp = polym_p[i] >= min_cov
        Sm = polym_m[i] >= min_cov
        
        for idx, j in enumerate(np.arange(0, 2*surr+1)-surr):
            SpS = _shift(Sp, j, False)
            SmS = _shift(Sm, -j, False)
            A_freq[idx] += np.sum(_shift(polym_p[i], j, 0)[A&SpS]) + np.sum(_shift(polym_m[i], -j, 0)[U&SmS])
            G_freq[idx] += np.sum(_shift(polym_p[i], j, 0)[G&SpS]) + np.sum(_shift(polym_m[i], -j, 0)[C&SmS])
            U_freq[idx] += np.sum(_shift(polym_p[i], j, 0)[U&SpS]) + np.sum(_shift(polym_m[i], -j, 0)[A&SmS])
            C_freq[idx] += np.sum(_shift(polym_p[i], j, 0)[C&SpS]) + np.sum(_shift(polym_m[i], -j, 0)[G&SmS])

    base_df = pd.DataFrame({'A': A_freq, 'U': U_freq, 'G': G_freq, 'C': C_freq})
    # normalize
    base_df = (base_df.T / base_df.T.sum()).T

    return(base_df)


def PRO_DMS_signal(filename, ref, start, end, strand, clip_filter=999):
    """
    Extract mutation and coverage matrices for CoSTseq analysis.

    Parameters:
    - filename (str): Path to the BAM file containing alignments.
    - ref (str): Reference sequence name.
    - start (int): Start position of the region.
    - end (int): End position of the region.
    - strand (str): Strand information ('+' or '-').
    - clip_filter (int): Maximum allowed number of clipped bases.

    Returns:
    - tuple: Mutation matrix and coverage matrix.
    """

    t_start = time.time()
    bamfile = pysam.AlignmentFile(filename, 'rb')

    mat_cov = np.zeros([end-start, end-start])
    mat_mut = np.zeros([end-start, end-start])

    try:        
        for read1, read2 in read_pair_generator(bamfile, ref, start+1, end+1):
            
            # strand check
            if (not read1.is_reverse and strand == '-') or (read1.is_reverse and strand == '+'):
                continue

            # determine RNA 3'-end position, which is always the 5'-end of read2
            _, mapped_5, mapped_3 = read_ends(read2)
            threeprime = mapped_5-start

            # 3'-end outside of requested region 
            if threeprime < 0 or end-mapped_5 <= 0: 
                continue
                
            # obtain mutational vectors, clipped bases
            read_map_1, clipped_3_1, clipped_5_1 = map_mutations(read1, with_clipped=True)
            read_map_2, clipped_3_2, clipped_5_2 = map_mutations(read2, with_clipped=True)

            # skip reads with 3'-end clipping >1 (5'-end of read2 is the 3'-end of the molecule)
            if len(clipped_5_2) > 1:
                continue

            # skip reads that have a specified amount of clipped bases elsewhere
            if (len(clipped_3_1) > clip_filter) or (len(clipped_5_1) > clip_filter) or (len(clipped_3_2) > clip_filter):
                continue 

            read_map = join_reads(read_map_1, read_map_2)

            # max. 10% of the read can be mutated
            if np.sum(np.isin(np.array(list(read_map.values())), ['A', 'G', 'C', 'T', '1'])) > 0.1 * len(read_map):
                continue

            # sort each nt into matrix
            for i, nt in enumerate(read_map.keys()):
                # skip if nt is outside requested region
                if nt-start < 0:
                    continue

                # skip overlap due to read2 trimming, stop if end of requested region is reached
                if (nt-start < threeprime and strand == '-') or (nt-start > threeprime and strand == '+') or (nt-start >= end-start):
                    continue

                # update mut and cov matrices
                if read_map[nt] in ['A', 'G', 'T', 'C', '1']:
                    mat_mut[threeprime, nt-start] += 1     
                
                if read_map[nt] in ['A', 'G', 'T', 'C', '1', '0']:
                    mat_cov[threeprime, nt-start] += 1

                # register the one allowed clipped base as mutation
                if strand == '+' and i == len(read_map)-1 and read_map[nt] == '.':
                    mat_mut[threeprime, nt-start] += 1
                    mat_cov[threeprime, nt-start] += 1
                elif strand == '-' and i == 0 and read_map[nt] == '.':
                    mat_mut[threeprime, nt-start] += 1
                    mat_cov[threeprime, nt-start] += 1

    finally:
        bamfile.close()
        t_end = time.time()
        print('time elapsed: %.2f min'%((t_end-t_start)/60))

    return(mat_mut, mat_cov)

def moving_average(x, w):
    """
    Calculate a simple moving average.

    Parameters:
    - x (numpy.ndarray): Input array.
    - w (int): Window size.

    Returns:
    - numpy.ndarray: Array with moving averages.
    """

    return (np.convolve(x, np.ones(w), 'valid') / w)


def binned_signal(signal_y, n_bins):
    """
    Bin a signal into specified number of bins.

    Parameters:
    - signal_y (numpy.ndarray): Input signal.
    - n_bins (int): Number of bins.

    Returns:
    - tuple: Normalized x-values, binned y-values, and bin size.
    """

    signal_x = np.arange(len(signal_y))/(len(signal_y)-1)
    norm_y, bin_edges, _ = binned_statistic(signal_x, signal_y, statistic='sum', bins=n_bins, range=None)
    bin_size = bin_edges[1] - bin_edges[0]
    norm_x = bin_edges[:-1]+0.5*bin_size
    
    return(norm_x, norm_y, bin_size)


def c_freq(seq):
    """
    Calculate the frequency of 'C' in a given RNA sequence.

    Parameters:
    - seq (str): RNA sequence.

    Returns:
    - numpy.ndarray: Frequency array.
    """
    
    freq = np.ones(len(seq))
    no_c = 0
    for i, nt in enumerate(seq):
        if nt.upper() == 'C':
            freq[i] += no_c
            no_c = 0
        else:
            no_c += 1
    return(freq)
