import pickle
import numpy as np
from src.bam_utils import reverse_complement
from copy import deepcopy
import pandas as pd
import matplotlib.colors as mplc
from matplotlib import cm
import matplotlib.pyplot as plt
from scipy.stats import binom
from tqdm import tqdm
import math
from scipy.stats import iqr

class CotrxMatrix():
    """
    Class representing cotranscriptional matrices.

    Methods:
    - __init__: Constructor method.
    - load_matrix_np: Load matrix from NumPy binary file.
    - load_matrix: Load matrix from pickle file.
    - write_matrix: Write matrix to pickle file.
    - process: Perform normalization, apply mask and remove positions that are not covered.
    - bin_matrix: Bin the matrix.
    - save_to_img: Save matrix to an image file.
    - get_residue_txn: Get reactivity trajectory for a specific nucleotide position.
    - calc_dynamic_pvals: Calculate dynamic p-values.
    - plot_avg_profile: Plot reactivity profile averaged over polymerase positions.
    - plot_single_pol_pos: Plot reactivity profile for single polymerase positions.
    """

    def __init__(self, filename, sample_id, seq, strand):
        """
        Constructor method for CotrxMatrix.

        Parameters:
        - filename (str): Path to the matrix file.
        - sample_id (str): Identifier for the sample.
        - seq (str): Sequence associated with the matrix.
        - strand (str): Strand information, either '+' or '-'.
        """

        self.sample_id = sample_id
        self.load_matrix(filename, strand)
        self.seq = seq
        
        self.mask_seq = np.array([True if nt in ['A', 'C'] else False for nt in self.seq])

    def load_matrix_np(self, filename, strand):
        """
        Load matrix from NumPy binary file. This is for compatibility and it's no longer used.

        Parameters:
        - filename (str): Path to the NumPy binary file.
        - strand (str): Strand information, either '+' or '-'.
        """
        with open(filename, 'rb') as f:
            self.mut, self.cov = pickle.load(f)
        self.size = self.mut.shape
        
        if strand == '-':
            self.mut = np.flip(self.mut, (0, 1))
            self.cov = np.flip(self.cov, (0, 1))

    def load_matrix(self, filename, strand):
        """
        Load matrix from a pickle file.

        Parameters:
        - filename (str): Path to the pickle file.
        - strand (str): Strand information, either '+' or '-'.
        """

        with open(filename, 'rb') as f:
            mut_d, cov_d, s = pickle.load(f)

        self.mut = dict_to_sparse(mut_d, s)
        self.cov = dict_to_sparse(cov_d, s)
        self.size = s

        if strand == '-':
            self.mut = np.flip(self.mut, (0, 1))
            self.cov = np.flip(self.cov, (0, 1))

    def write_matrix(self, filename):
        """
        Write matrix to a pickle file.

        Parameters:
        - filename (str): Output path for the pickle file.
        """
        mut_d = sparse_to_dict(self.mut)
        cov_d = sparse_to_dict(self.cov)

        with open(filename, 'wb') as f:
            pickle.dump((mut_d, cov_d, self.size), f)

    def process(self, norm=True, reassign=False, max_dist=15, per_exclude=10, min_cov=800):
        """
        Process the matrix.

        Parameters:
        - norm (bool): Whether to normalize the matrix (default is True).
        - reassign (bool): Whether to reassign polymerase positions to nearest run-on nucleotide (default is False).
        - max_dist (int): Maximum distance for reassignment (default is 15).
        - per_exclude (int): Upper percentile of values to exclude during normalization (default is 10).
        - min_cov (int): Minimum coverage threshold (default is 800).
        """

        if reassign:
            self.mut, self.cov = reassign_pols(self.mut, self.cov, self.seq, max_dist=max_dist)

        if norm:
            self.rea, self.fac_norm = normalize(process_folding_matrix(self.mut, self.cov, min_cov=min_cov), self.mask_seq, per_exclude=per_exclude)
            self.rea_full, _ = normalize(process_folding_matrix(self.mut, self.cov, return_full=True, min_cov=min_cov), self.mask_seq, per_exclude=per_exclude)
        else:
            self.rea = process_folding_matrix(self.mut, self.cov, min_cov=min_cov)
            self.rea_full = process_folding_matrix(self.mut, self.cov, return_full=True, min_cov=min_cov)
            self.fac_norm = 1

    def bin_matrix(self, binsize, **kwargs):
        """
        Bin the matrix.

        Parameters:
        - binsize (int): Size of the bins for binning the matrix.
        - **kwargs: Additional keyword arguments for processing the folding matrix.
        """

        # reactivity calculation in discrete bins of polymerase positions
        bins = np.arange(0, self.size[0], binsize)
        mut_bin = np.zeros([len(bins), self.size[0]])
        cov_bin = np.zeros([len(bins), self.size[0]])

        for i in range(len(bins)):
            if i < len(bins)-1:
                cov_bin[i, :] = np.sum(self.cov[bins[i]:bins[i+1], :], axis=0)
                mut_bin[i, :] = np.sum(self.mut[bins[i]:bins[i+1], :], axis=0)

            else: # last bin in case it doesn't add up
                cov_bin[i, :] = np.sum(self.cov[bins[i]:, :], axis=0)
                mut_bin[i, :] = np.sum(self.mut[bins[i]:, :], axis=0)
                
        self.bin_mut = mut_bin
        self.bin_cov = cov_bin
        self.bin_rea = process_folding_matrix(mut_bin, cov_bin, **kwargs)


    def save_to_img(self, out_name, mask=None, filter_ac=True, **kwargs):
        """
        Save matrix to an image file.

        Parameters:
        - out_name (str): Output file name.
        - mask (numpy.ndarray): Mask for filtering values (default is None).
        - filter_ac (bool): Whether to show only A and C nucleotides (default is True).
        - **kwargs: Additional keyword arguments for saving the image.
        """

        if filter_ac and mask is None:
            save_mat_to_img(self.rea, out_name, mask=self.mask_seq, **kwargs)
        elif mask is not None:
            save_mat_to_img(self.rea, out_name, mask=mask, **kwargs)
        else:
            save_mat_to_img(self.rea, out_name, **kwargs)


    def get_residue_txn(self, res, min_cov=800, binned=False):
        """
        Get reactivity for a specific nucleotide position.

        Parameters:
        - res (int): Nucleotide position.
        - min_cov (int): Minimum coverage threshold (default is 800).
        - binned (bool): Whether to use binned data (default is False).

        Returns:
        - Tuple[numpy.ndarray, numpy.ndarray]: Tuple containing polymerase positions and reactivities.
        """
        
        if binned:
            y = process_folding_matrix(self.bin_mut, self.bin_cov, return_full=True, min_cov=min_cov, min_frac=0)[:, res]
        else:
            y = self.rea_full[:, res]
        
        x = np.argwhere(~np.isnan(y)) + 1
        y = y[~np.isnan(y)]
        
        return(x.squeeze(), y)

    def calc_dynamic_pvals(self, min_N=10, min_cov=800):
        """
        Calculate dynamic p-values. This never really worked.

        Parameters:
        - min_N (int): Minimum number of data points for p-value calculation (default is 10).
        - min_cov (int): Minimum coverage threshold (default is 800).
        """
    
        cov = self.cov.astype(float)
        mut = self.mut.astype(float)
        
        mat_mask = cov >= min_cov
        filt_cov = np.sum(mat_mask, axis=1)
        
        cov[~mat_mask] = np.nan
        mut[~mat_mask] = np.nan
        
        self.pvals = []
        
        for nt in range(self.size[0]):

            y_mut = mut[:, nt]
            y_cov = cov[:, nt]

            y_mut = y_mut[~np.isnan(y_mut)].astype(int)
            y_cov = y_cov[~np.isnan(y_cov)].astype(int)

            if len(y_mut) >= min_N:
                self.pvals.append(calc_p_rea_change(y_mut, y_cov))
            else:
                self.pvals.append(np.nan)
        
        self.pvals = np.array(self.pvals)

    def plot_avg_profile(self, start, stop, cmap_loc='cmap.txt'):
        """
        Plot average reactivity profile.

        Parameters:
        - start (int): Starting position for the plot.
        - stop (int): Ending position for the plot.
        - cmap_loc (str): Location of the colormap file (default is 'cmap.txt').
        """

        cmap = np.array(pd.read_csv(cmap_loc, header=None))/255
        cmap = mplc.ListedColormap(cmap)
        norm = mplc.Normalize(vmin=0, vmax=1)
        m = cm.ScalarMappable(norm=norm, cmap=cmap)

        fig = plt.figure(figsize=(18, 3))

        x = np.arange(start, stop)
        y = np.nanmean(self.rea, axis=0)[start:stop]

        clrs = m.to_rgba(y)
        clrs[~self.mask_seq[start:stop]] = [230/255, 230/255, 230/255, 1]
        
        plt.bar(x, y, 1, color=clrs)
        plt.colorbar(m)
        plt.ylim([0, 3])
        plt.show()

    def plot_single_pol_pos(self, pos, min_cov=800, cmap_loc='cmap.txt'):
        """
        Plot reactivity for single polymerase positions.

        Parameters:
        - pos (List[int]): List of polymerase positions to plot.
        - min_cov (int): Minimum coverage threshold (default is 800).
        - cmap_loc (str): Location of the colormap file (default is 'cmap.txt').

        Returns:
        - Dict[int, Tuple[numpy.ndarray, numpy.ndarray]]: Dictionary containing positions and reactivities.
        """

        cmap = np.array(pd.read_csv(cmap_loc, header=None))/255
        cmap = mplc.ListedColormap(cmap)
        norm = mplc.Normalize(vmin=0, vmax=1)
        m = cm.ScalarMappable(norm=norm, cmap=cmap)

        fig, axs = plt.subplots(len(pos), 1, figsize=(18, 1*len(pos)), sharex=True)

        rea_dict = {}
        for i, p in enumerate(pos):

            y = self.rea_full[p,:]
            x = np.arange(1, len(y)+1)

            clrs = m.to_rgba(y)
            clrs[self.cov[p,:] < min_cov] = [1, 0, 0, 1]
            clrs[~self.mask_seq] = [230/255, 230/255, 230/255, 1]
            
            axs[i].bar(x, y, 1, color=clrs)
            rea_dict[p] = (x, y)

        plt.setp(axs, ylim=[0, 3])
        plt.show()
        
        return(rea_dict)



def combine_reps(matrices, sample_id):
    """
    Combine multiple CotrxMatrix objects into a single matrix.

    Parameters:
    - matrices (List[CotrxMatrix]): List of CotrxMatrix objects to be combined.
    - sample_id (str): Sample identifier for the new matrix.

    Returns:
    - CotrxMatrix: Combined CotrxMatrix object.
    """

    assert isinstance(matrices, list), "Input must be list of CotrxMatrix objects."
    
    new = deepcopy(matrices[0])
    if hasattr(new, 'rea'):
        delattr(new, 'rea')
    new.cov = np.sum([matrix.cov for matrix in matrices], axis=0)
    new.mut = np.sum([matrix.mut for matrix in matrices], axis=0)

    new.sample_id = sample_id
    return(new)

def reassign_pols(mut, cov, seq, max_dist=5, nts_runon=['C']):
    """
    Reassign reactivity vectors to the nearest run-on nucleotide.

    Parameters:
    - mut (numpy.ndarray): Matrix of mutations.
    - cov (numpy.ndarray): Matrix of coverage.
    - seq (str): Nucleotide sequence.
    - max_dist (int): Maximum distance to search for valid C positions (default is 5).
    - nts_runon (List[str]): List of nucleotides that indicate run-on nucleotides (default is ['C']).

    Returns:
    - Tuple[numpy.ndarray, numpy.ndarray]: Mutations and coverage matrices after reassignment.
    """

    def _find_next_valid(i, seq, max_dist, nts_runon):
        for d in range(1, max_dist+1):
            if seq[i-d] in nts_runon:
                return(i-d)

    mut_re = np.zeros(mut.shape)
    cov_re = np.zeros(cov.shape)
    nts_lost = 0
    for i, nt in enumerate(seq):
        if nt in nts_runon:
            mut_re[i,:] += mut[i,:]
            cov_re[i,:] += cov[i,:]
        else:
            # reassign pol position to previous valid nucleotide
            j = _find_next_valid(i, seq, max_dist, nts_runon)
            if j is not None:
                mut_re[j,:] += mut[i,:]
                cov_re[j,:] += cov[i,:]
            else:
                nts_lost += 1

    print(f"{nts_lost} Pol positions lost during reassignment.")
    return(mut_re, cov_re)


def save_mat_to_img(mat, out_name, coords=None, mask=None, cmap_loc='cmap.txt', vmin=0, vmax=1):
    """
    Save matrix to an image file.

    Parameters:
    - mat (numpy.ndarray): Matrix to be saved.
    - out_name (str): Output file name.
    - coords (Tuple[int, int]): Coordinates of the matrix to be saved (default is None).
    - mask (numpy.ndarray): Mask to filter the matrix (default is None).
    - cmap_loc (str): Location of the colormap file (default is 'cmap.txt').
    - vmin (float): Minimum value for colormap normalization (default is 0).
    - vmax (float): Maximum value for colormap normalization (default is 1).
    """
        
    cmap = np.array(pd.read_csv(cmap_loc, header=None))/255
    cmap = mplc.ListedColormap(cmap)
    c_norm = plt.Normalize(vmin=vmin, vmax=vmax)
    
    if coords is not None:
        mat = mat[:, coords[0]:coords[1]]
        if mask is not None:
            mask = mask[coords[0]:coords[1]]

    image = cmap(c_norm(mat))
    if mask is not None:
        image[:, ~mask] = [230/255, 230/255, 230/255, 1]

    plt.imsave(out_name, image)
    plt.imshow(image)


def sparse_to_dict(mat):
    """
    Convert a sparse matrix to a dictionary.

    Parameters:
    - mat (numpy.ndarray): Sparse matrix.

    Returns:
    - Dict[Tuple[int, int], int]: Dictionary representation of the sparse matrix.
    """

    mat_d = {}
    s = mat.shape
    for i in range(s[0]):
        for j in range(s[1]):
            if mat[i,j] > 0:
                mat_d[i, j] = mat[i,j]
    return(mat_d)

def dict_to_sparse(mat_d, s, of_type=int):
    """
    Convert a dictionary to a sparse matrix.

    Parameters:
    - mat_d (Dict[Tuple[int, int], int]): Dictionary representation of the matrix.
    - s (Tuple[int, int]): Size of the matrix.
    - of_type (Type): Data type of the matrix (default is int).

    Returns:
    - numpy.ndarray: Sparse matrix.
    """

    mat = np.zeros(s, dtype=of_type)
    for i, j in mat_d.keys():
        mat[i, j] = mat_d[i, j]
    return(mat)

def process_folding_matrix(mut, cov, return_full=False, min_cov=800, min_frac=0.01):
    """
    Process folding matrix.

    Parameters:
    - mut (numpy.ndarray): Matrix of mutations.
    - cov (numpy.ndarray): Matrix of coverage.
    - return_full (bool): Whether to return the full matrix or to remove empty rows (default is False).
    - min_cov (int): Minimum coverage threshold (default is 800).
    - min_frac (float): Minimum fraction of zero coverage for filtering (default is 0.01).

    Returns:
    - numpy.ndarray: 
    """

    # coverage filter
    mat_mask = cov >= min_cov
    filt_cov = np.sum(mat_mask, axis=1) >= min_frac*np.sum(cov==0, axis=1)
  
    if not return_full:
        rea = (mut/cov)[filt_cov]
        rea[~mat_mask[filt_cov]] = np.nan
    else:
        rea = mut/cov
        rea[~mat_mask] = np.nan

    return(rea)

def calc_p_rea_change(vec_mut, vec_cov):
    """
    Calculate p-value for reactivity change.

    Parameters:
    - vec_mut (numpy.ndarray): Vector of mutations.
    - vec_cov (numpy.ndarray): Vector of coverage.

    Returns:
    - float: Calculated p-value for reactivity change.
    """
    return(np.prod( binom.pmf( vec_mut, vec_cov, np.mean(vec_mut/vec_cov) ) ) ** ( 1/len(vec_mut) ))


def normalize(arr, mask=None, per_exclude=10):
    """
    Normalize reactivities in co-transcriptional folding matrix.

    Parameters:
    - arr (numpy.ndarray): Array to be normalized.
    - mask (numpy.ndarray): Mask for filtering values (default is None).
    - per_exclude (int): Upper percentile of values to exclude during normalization (default is 10).

    Returns:
    - Tuple[numpy.ndarray, float]: Normalized array and normalization factor.
    """

    # sort and remove nan
    if mask is not None:
        if len(arr.shape) > 1:
            arr_sort = arr[:, mask]
        else:
            arr_sort = arr[mask]
        arr_sort = np.sort(arr_sort[arr_sort==arr_sort]).flatten()
    else:
        arr_sort = np.sort(arr[arr==arr]).flatten()

    # define values excluded for normalization
    IQR = iqr(arr_sort)
    PER = arr_sort[arr_sort.shape[0] - 1 - (len(arr_sort)//per_exclude if len(arr_sort) > 100 else len(arr_sort)//5)]
    threshold = np.max([1.5*IQR, PER])
    good = (arr_sort < threshold)
    
    # calculate normalization factor (mean of top 10 percent -> length according to ORIGINAL vector)
    # I think length should be according to the new vector (after excluding values), so potentially change this in the future
    top_rea = arr_sort[good][(-len(arr_sort)//10)+1:]
    fac_norm = 1/np.mean(top_rea)
    
    return(arr*fac_norm, fac_norm)


def write_HDProbe(out_name, matrices, labels, seq, binned=False):
    """
    Write co-transcriptional matrices to HDProbe format file.

    Parameters:
    - out_name (str): Output file name.
    - matrices (List[CotrxMatrix]): List of CotrxMatrix objects.
    - labels (List[str]): List of labels for each matrix.
    - seq (str): Sequence.
    - binned (bool): Whether the matrices are binned (default is False).
    """

    assert isinstance(matrices, list), "Input must be list of CotrxMatrix objects."
    
    with open(out_name, 'w') as f:
        # header
        f.write('\t'.join(['sample_name', 'chromosome', 'genomic_position', 'nt', 'strand', 'n_mutations', 'n_reads' + '\n']))
        
        for matrix, label in tqdm(zip(matrices, labels), total=len(matrices)):
            if not binned:
                cov_d = sparse_to_dict(matrix.cov)
                mut_d = sparse_to_dict(matrix.mut)
            else:
                cov_d = sparse_to_dict(matrix.bin_cov)
                mut_d = sparse_to_dict(matrix.bin_mut)
            
            for pos in mut_d.keys():
                f.write('\t'.join([label, str(pos[0]), str(pos[1]), seq[pos[1]], '+', str(mut_d[pos]), str(cov_d[pos]) + '\n']))

def load_HDProbe(filename, s):
    """
    Load HDProbe format file.

    Parameters:
    - filename (str): Input file name.
    - s (Tuple[int, int]): Size of the matrix.

    Returns:
    - numpy.ndarray: Loaded matrix.
    """

    hdp = pd.read_csv(filename)[['chr', 'pos', 'padj_BH']]
    pval_d = {}
    for i, nt in hdp.iterrows():
        if nt['chr'] == nt['chr'] and nt['pos'] == nt['pos'] and nt['padj_BH'] == nt['padj_BH']: # escape nts that aren't matched between data sets
            pval_d[int(nt['chr']), int(nt['pos'])] = nt['padj_BH']
        else:
            print('Found NaN in data.')

    pval_mat = dict_to_sparse(pval_d, s, of_type=float)
    pval_mat[pval_mat==0] = np.nan

    return(pval_mat)