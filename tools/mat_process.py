""""From https://github.com/cmdoret/hicreppy"""

from typing import Union

import numpy as np
import scipy.sparse as sp


def smooth(signal: "sp.coo_matrix", h: int) -> "sp.coo_matrix":
    """
    Smooth (blur) input sparse Hi-C sparse matrix using a uniform kernel of width 2*h+1.

    Parameters
    ----------
    signal : scipy.sparse.coo_matrix
        Input intrachromosomal Hi-C matrix to smooth (upper triangle matrix).
    h : int
        Half width of the kernel (neighbourhood size).

    Returns
    -------
    scipy.sparse.coo_matrix :
        The smoothed matrix.
    """
    # km = 2 * h + 1
    # kernel = np.ones((km, km)) / (km**2)
    # sm, sn = signal.shape

    # # Sanity checks
    # if sp.issparse(kernel):
    #     raise ValueError("cannot handle kernel in sparse format")
    # if not sp.issparse(signal):
    #     raise ValueError("cannot handle signal in dense format")
    # constant_kernel = kernel[0, 0]

    # # Simplified convolution for the special case where kernel is constant:
    # l_subkernel_sp = sp.diags(np.ones(km), np.arange(km), shape=(sm - km + 1, sm), format="csr")
    # r_subkernel_sp = sp.diags(np.ones(km), -np.arange(km), shape=(sn, sn - km + 1), format="csr")
    # out = (l_subkernel_sp @ signal) @ r_subkernel_sp
    # out *= constant_kernel
    # # Resize matrix: increment rows and cols by half kernel and set shape to input
    # # matrix, effectively adding margins.
    # out = out.tocoo()
    # rows, cols = out.row + h, out.col + h
    # out = sp.coo_matrix((out.data, (rows, cols)), shape=(sm, sn), dtype=np.float64)

    # return out

    # pad signal matrix with h zeros all around
    rows, cols = signal.row + h, signal.col + h
    sm, sn = signal.shape
    signal = sp.coo_matrix((signal.data, (rows, cols)), shape=(sm + 2 * h, sn + 2 * h), dtype=np.float64)

    sm, sn = signal.shape

    km = 2 * h + 1
    kernel = np.ones((km, km)) / (km ** 2)
    sm, sn = signal.shape

    # Sanity checks
    if sp.issparse(kernel):
        raise ValueError("cannot handle kernel in sparse format")
    if not sp.issparse(signal):
        raise ValueError("cannot handle signal in dense format")
    constant_kernel = kernel[0, 0]

    # Simplified convolution for the special case where kernel is constant:
    l_subkernel_sp = sp.diags(np.ones(km), np.arange(km), shape=(sm - km + 1, sm), format="csr")
    r_subkernel_sp = sp.diags(np.ones(km), -np.arange(km), shape=(sn, sn - km + 1), format="csr")
    out = (l_subkernel_sp @ signal) @ r_subkernel_sp
    out *= constant_kernel
    # Resize matrix: increment rows and cols by half kernel and set shape to input
    # matrix, effectively adding margins.
    out = out.tocoo()

    return out


def vstrans(d1: "np.ndarray[float]", d2: "np.ndarray[float]") -> "np.ndarray[float]":
    """
    Variance stabilizing transformation.

    Variance stabilizing transformation to normalize read counts before
    computing stratum correlation. This normalizes counts so that different
    strata share similar dynamic ranges.

    Parameters
    ----------
    d1 : numpy.ndarray of floats
        Diagonal of the first matrix.
    d2 : numpy.ndarray of floats
        Diagonal of the second matrix.

    Returns
    -------
    r2k : numpy.ndarray of floats
        Array of weights to use to normalize counts.
    """
    # Get ranks of counts in diagonal
    ranks_1 = np.argsort(d1) + 1
    ranks_2 = np.argsort(d2) + 1
    # Scale ranks betweeen 0 and 1
    nranks_1 = ranks_1 / max(ranks_1)
    nranks_2 = ranks_2 / max(ranks_2)
    nk = len(ranks_1)
    r2k = np.sqrt(np.var(nranks_1 / nk) * np.var(nranks_2 / nk))
    return r2k


def subsample_contacts(M: "sp.coo_matrix", n_contacts: float) -> "sp.coo_matrix":
    """Bootstrap sampling of contacts in a sparse Hi-C map.

    Parameters
    ----------
    M : scipy.sparse.coo_matrix
        The input Hi-C contact map in sparse format.
    n_contacts : float
        The number of contacts to be sampled if larger than one.
        The proportion of contacts to be sampled if between 0 and 1.

    Returns
    -------
    scipy.sparse.coo_matrix
        A new matrix with a fraction of the original contacts.
    """
    try:
        if n_contacts < 0:
            raise ValueError("n_contacts must be strictly positive")
        if n_contacts <= 1:
            n_contacts *= M.data.sum()
    except ValueError:
        raise ValueError("n_contacts must be a float")
    S = M.data.copy()
    # Match cell idx to cumulative number of contacts
    cum_counts = np.cumsum(S)
    # Total number of contacts to sample
    tot_contacts = int(cum_counts[-1])

    # Sample desired number of contacts from the range(0, n_contacts) array
    sampled_contacts = np.random.choice(int(tot_contacts), size=int(n_contacts), replace=False)

    # Get indices of sampled contacts in the cum_counts array
    idx = np.searchsorted(cum_counts, sampled_contacts, side="right")

    # Bin those indices to the same dimensions as matrix data to get counts
    sampled_counts = np.bincount(idx, minlength=S.shape[0])

    # Get nonzero values to build new sparse matrix
    nnz_mask = sampled_counts > 0
    sampled_counts = sampled_counts[nnz_mask].astype(np.float32)
    sampled_rows = M.row[nnz_mask]
    sampled_cols = M.col[nnz_mask]

    return sp.coo_matrix((sampled_counts, (sampled_rows, sampled_cols)), shape=(M.shape[0], M.shape[1]),)


def diag_trim(mat: Union["sp.dia_matrix", "np.ndarray"], n: int) -> Union["sp.dia_matrix", "np.ndarray"]:
    """
    Trim an upper triangle sparse matrix so that only the first n diagonals are kept.

    Parameters
    ----------
    mat : sp.dia_matrix or numpy.ndarray
        The sparse matrix to be trimmed
    n : int
        The number of diagonals from the center to keep (0-based).

    Returns
    -------
    sp.dia_matrix or numpy.ndarray:
        The diagonally trimmed upper triangle matrix with only the first n
        diagonal.
    """
    if not sp.issparse(mat):
        trimmed = mat.copy()
        n_diags = trimmed.shape[0]
        for diag in range(n, n_diags):
            set_mat_diag(trimmed, diag, 0)
        return trimmed
    if mat.format != "dia":
        raise ValueError("input type must be sp.dia_matrix")
    # Create a new matrix from the diagonals below max dist (faster than removing them)
    keep_offsets = np.flatnonzero((mat.offsets <= n) & (mat.offsets >= 0))

    trimmed = sp.dia_matrix((mat.data[keep_offsets], mat.offsets[keep_offsets]), shape=mat.shape)

    return trimmed


def set_mat_diag(mat: "np.ndarray", diag: int = 0, val: int = 0) -> float:
    """
    Set the nth diagonal of a symmetric 2D numpy array to a fixed value.
    Operates in place.

    Parameters
    ----------
    mat : numpy.ndarray
        Symmetric 2D array of floats.
    diag : int
        0-based index of the diagonal to modify. Use negative values for the
        lower half.
    val : float
        Value to use for filling the diagonal
    """
    m = mat.shape[0]
    start = diag
    end = m ** 2 - diag * m
    step = m + 1
    mat.flat[start:end:step] = val
