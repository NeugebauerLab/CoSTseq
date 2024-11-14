import numpy as np
import matplotlib.pyplot as plt

def dscatter_plot(X, Y, nbins=[], order=False, lamb=20, markersize=5, ax=None, **kwargs):
    """Scatter plot with smoothed densities.

    Arguments:
    X, Y -- data
    nbins -- number of bins in X and Y direction. Tuple-like.
    order -- Should the highest-density points be plotted on tip? Decreases performance.
    lamb -- smoothness parameter
    markersize -- size of points
    ax -- axis to plot on, e.g. if using plt.subplots
    kwargs -- other plt.scatter() arguments
    """

    col, _, _, _ = pydscatter(X, Y, nbins, lamb)

    if order:
        sorting = np.argsort(col)
    else:
        sorting = np.arange(0, len(col))

    if not ax:        
        p = plt.scatter(np.array(X)[sorting], np.array(Y)[sorting], markersize, col[sorting], **kwargs)
        return(p)
    else:
        p = ax.scatter(np.array(X)[sorting], np.array(Y)[sorting], markersize, col[sorting], **kwargs)
        return(p)

def dscatter_img(X, Y, nbins=[], lamb=20, ax=None, **kwargs):
    """Rendered image of smoothes densities.

    Arguments:
    X, Y -- data
    nbins -- number of bins in X and Y direction. Tuple-like.
    lamb -- smoothness parameter
    ax -- axis to plot on, e.g. if using plt.subplots
    kwargs -- other plt.imshow() arguments
    """
    
    _, F, ctrs1, ctrs2 = pydscatter(X, Y, nbins, lamb)

    if not ax:        
        p = plt.imshow(np.flipud(F), extent=(min(ctrs1), max(ctrs1), min(ctrs2), max(ctrs2)), aspect='auto', **kwargs)
        return(p)
    else:
        p = ax.imshow(np.flipud(F), extent=(min(ctrs1), max(ctrs1), min(ctrs2), max(ctrs2)), aspect='auto', **kwargs)
        return(p)

def dscatter_contour(X, Y, nbins=[], lamb=20, ax=None, **kwargs):
    """Contour plot of densities.

    Arguments:
    X, Y -- data
    nbins -- number of bins in X and Y direction. Tuple-like.
    lamb -- smoothness parameter
    ax -- axis to plot on, e.g. if using plt.subplots
    kwargs -- other plt.contour() arguments
    """

    col, F, cntr1, cntr2 = pydscatter(X, Y, nbins, lamb)

    if not ax:        
        p = plt.contour(cntr1, cntr2, F, **kwargs)
        return(p)
    else:
        p = ax.contour(cntr1, cntr2, F, **kwargs)
        return(p)

def pydscatter(X, Y, nbins=[], lamb=20):
    """Calculated density values based on data.

    Arguments:
    X, Y -- data
    nbins -- number of bins in X and Y direction. Tuple-like.
    lamb -- smoothness parameter
    
    Returns:
    col -- density values
    F -- density map
    ctrs1, ctrs2 -- contour coordinates for contour plot
    """

    if not nbins:
        nbins = [min(len(np.unique(X)), 200), min(len(np.unique(Y)), 200)]

    edges1 = np.linspace(min(X), max(X), nbins[0]+1)
    ctrs1 = edges1[:-1] + 0.5*np.diff(edges1)
    edges2 = np.linspace(min(Y), max(Y), nbins[1]+1);
    ctrs2 = edges2[:-1] + 0.5*np.diff(edges2)

    bins = np.vstack((np.digitize(Y, edges2), np.digitize(X, edges1)))

    H = aggregate(bins, 1, size=nbins[::-1])
    G = smooth1D(H, nbins[1]/lamb)
    F = smooth1D(G.conjugate().T, nbins[0]/lamb).conjugate().T

    I = F/np.max(F)
    ind = sub2ind(I.shape, bins[0,:], bins[1,:])
    col = I.flatten()[ind]

    return(col, F, ctrs1, ctrs2)

def aggregate(group_idx, a, size):
    """Similar to MATLAB accumarray"""

    group_idx = np.ravel_multi_index(group_idx, size, mode='clip')    
    ret = np.bincount(group_idx, minlength=size[0] * size[1])
    if a != 1:
        ret *= a
    ret = ret.reshape(size)
    return(ret)

def smooth1D(Y, lamb):
    m = Y.shape[0]
    E = np.eye(m)
    D1 = np.diff(E)
    D2 = np.diff(D1)
    P = lamb**2 * np.dot(D2, D2.conjugate().T) + 2*lamb * np.dot(D1, D1.conjugate().T)
    Z,_,_,_ = np.linalg.lstsq(E+P, Y, rcond=None)
    return(Z)

def sub2ind(array_shape, rows, cols):
    """Same as MATLAB sub2int"""
    ind = rows*array_shape[1] + cols
    ind[ind < 0] = -1
    ind[ind >= array_shape[0]*array_shape[1]] = -1
    return(ind)