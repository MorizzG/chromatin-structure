#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 15:47:59 2022.

@author: mg

SSC calculation based on https://github.com/cmdoret/hicreppy
"""
import warnings
from multiprocessing import Pool
from typing import Optional

import numpy as np
import scipy.sparse as sp
import scipy.stats as ss

from . import mat_process


def _sample_scc(args):
    chrom_1, chrom_2, h_value, max_bins = args
    sub_1 = mat_process.subsample_contacts(chrom_1, 0.1)
    sub_2 = mat_process.subsample_contacts(chrom_2, 0.1)
    # Smooth the matrix using a mean filter
    smooth_1 = mat_process.smooth(sub_1, h_value)
    smooth_2 = mat_process.smooth(sub_2, h_value)
    return get_scc(smooth_1, smooth_2, max_bins=max_bins)


def h_train(
    mat1: "sp.csr_matrix",
    mat2: "sp.csr_matrix",
    binsize: int,
    chroms_lengths: "np.array[int]",
    max_dist: int,
    h_max: int,
    sample_count: int = 100,
    full_range: bool = False,
) -> int:
    """
    Find the optimal value for the smoothing parameter h.

    For each value of h
    in h_range, in ascending order, separate intrachromosomal sub-matrices of
    each chromosome. Each sub-matrix is subsampled to 10% contacts and the
    stratum adjusted correlation coefficient (SCC) is computed 10 times between
    the corresponding chromosome of both samples. The mean of those 10
    repetitions is computed for each chromosome and their mean weighted by
    chromosome length is used. The optimal h value is defined as the smallest
    value for which the increment of scc is less than 0.01.


    Parameters
    ----------
    mat1 : cooler.Cooler
        First matrix to compare.
    mat2 : cooler.Cooler
        Second matrix to compare.
    max_dist : int
        Maximum distance at which ton consider interactions, in basepairs.
    h_max : int
        The maximum value of the smoothing parameter h (neighbourhood size) to
        test. All values between 0 and h_max will be explored.
    sample_count : int
        Number of samples to draw for each h-value

    Returns
    -------
    int :
        Optimal value for the smoothing parameter h.
    """
    print("Starting h-training\n")
    max_bins = max_dist // binsize
    # Define chromosomes to scan
    # NOTE: chromosomes smaller than the kernel used for smoothing or the max_dist must
    # be removed.
    # min_size = max((2 * h_max + 1) * mat1.binsize, max_dist)
    # chromlist, chroms_lengths = make_chromlist(
    #     mat1, whitelist, blacklist, min_size=min_size
    # )
    # if not len(chromlist):
    #     raise KeyError(
    #         "All chromosomes were too short and have been discarded. "
    #         "Try reducing max_dist."
    #     )
    chrom_starts = np.insert(np.cumsum(chroms_lengths), 0, 0)
    prev_scc = -np.inf
    sccs = []
    p = Pool(8)
    for h_value in range(h_max + 1):
        # Compute SCC values separately for each chromosome
        chroms_scc = np.zeros(len(chroms_lengths))
        for c, chrom in enumerate(chroms_lengths):
            # samples_scc = np.zeros(10)
            c_start = chrom_starts[c]
            c_end = chrom_starts[c + 1]
            chrom_1 = mat1[c_start:c_end, c_start:c_end]
            chrom_2 = mat2[c_start:c_end, c_start:c_end]
            # Trim diagonals which are too far to be scanned to reduce
            # compute time and memory usage
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=sp.SparseEfficiencyWarning)
                chrom_1 = mat_process.diag_trim(chrom_1.todia(), max_bins + h_value).tocoo()
                chrom_2 = mat_process.diag_trim(chrom_2.todia(), max_bins + h_value).tocoo()
            # Sample 10% contacts and smooth 10 times for this chromosome
            samples_scc = np.array(p.map(_sample_scc, sample_count * [(chrom_1, chrom_2, h_value, max_bins)]))
            # Use average SCC from 10 subsamples
            chroms_scc[c] = np.nanmean(samples_scc)
        # Compute the genome SCC for this value of h using the weighted averge
        # of chromosomes SCC by their lengths. NaN values of SCC are not considered
        # This happens when comparing empty diagonals
        nan_scc_mask = ~np.isnan(chroms_scc)
        trunc_scc = chroms_scc[nan_scc_mask]
        trunc_lengths = np.array(chroms_lengths)[nan_scc_mask]
        curr_scc = np.average(trunc_scc, weights=trunc_lengths)
        sccs += [curr_scc]
        print(f"Found SCC of {curr_scc:.3f} with h={h_value:>2d} " f"with improvement {curr_scc - prev_scc:.3f}")
        # Check if SCC improvement is less than threshold
        delta = curr_scc - prev_scc
        if not full_range and 0 < delta and delta < 0.01:
            break
        prev_scc = curr_scc
    if full_range:
        return h_value, sccs

    if h_value == h_max:
        print("Note: It's likely that your searching range is too " "narrow. Try to expand the range and rerun it.\n",)
        return h_value, sccs
    # Return last h value with improvement >= 0.01
    # print(f"{sccs=}")
    return max(h_value - 1, 0), sccs


def get_scc(mat1: "sp.csr_matrix", mat2: "sp.csr_matrix", max_bins: int) -> float:
    """
    Calculate SCC of two matrices.

    Compute the stratum-adjusted correlation coefficient (SCC) between two
    Hi-C matrices up to max_dist. A Pearson correlation coefficient is computed
    for each diagonal in the range of 0 to max_dist and a weighted sum of those
    coefficients is returned.

    Parameters
    ----------
    mat1 : sp.csr_matrix
        First matrix to compare.
    mat2 : sp.csr_matrix
        Second matrix to compare.
    max_bins : int
        Maximum distance at which to consider, in bins.

    Returns
    -------
    scc : float
        Stratum adjusted correlation coefficient.
    """
    corr_diag = np.zeros(len(range(max_bins)))
    weight_diag = corr_diag.copy()
    for d in range(max_bins):
        d1 = mat1.diagonal(d)
        d2 = mat2.diagonal(d)
        # Silence NaN warnings: this happens for equal diagonals
        # which will be dropped
        with warnings.catch_warnings():
            # Warning does not exist in older scipy versions (<=1.2)
            try:
                warnings.filterwarnings("ignore", category=ss.PearsonRConstantInputWarning)
            except AttributeError:
                pass
            # Compute raw pearson coeff for this diag
            r = ss.pearsonr(d1, d2)[0]
            corr_diag[d] = r if not np.isnan(r) else 1
        # Compute weight for this diag
        r2k = mat_process.vstrans(d1, d2)
        weight_diag[d] = len(d1) * r2k
    # Normalize weights
    weight_diag /= sum(weight_diag)

    # Weighted sum of coefficients to get SCCs
    # scc = np.sum(corr_diag * weight_diag)
    scc = np.average(corr_diag, weights=weight_diag)

    return scc


def genome_scc(
    mat1: "sp.csr_matrix",
    mat2: "sp.csr_matrix",
    binsize: int,
    chroms_lengths: "np.array[int]",
    max_dist: int,
    h: int,
    subsample: Optional[float] = None,
) -> float:
    """Compute the Stratum-adjusted correlation coefficient (SCC) for the whole genome.

    Separate intrachromosomal sub-matrices of each chromosome. Compute the
    stratum adjusted correlation coefficient (SCC) between the corresponding
    chromosome of both samples.  The mean of SCCs weighted by chromosome length
    is used.

    Parameters
    ----------
    mat1 : sp.csr_matrix
        First matrix to compare.
    mat2 : sp.csr_matrix
        Second matrix to compare.
    binsize : int
        binsize of the matrices
    chroms_lengths : pandas.Series[int]
        list of chromosome lengths
    max_dist : int
        Maximum distance at which ton consider interactions, in basepairs.
    h : int
        The smoothing parameter. A mean filter is used to smooth matrices, this
        value is the size of the filter.
    subsample : None or float
        The number of contacts to which matrices should be subsampled, if
        greater than 1; the proportion of contacts to keep, if between 0 and 1.
        When you plan to compare multiple matrices, it can be useful to
        subsample all of them to the same value to remove potential biases
        caused by different coverages. Set to 0 to disable subsampling.

    Returns
    -------
    scc : float
        Stratum adjusted correlation coefficient.
    """
    if type(mat1) is not sp.csr_matrix or type(mat2) is not sp.csr_matrix:
        raise ValueError("Matrices must be scipy.sparse CSR matrices")
    # if subsample is not None:
    #     raise NotImplementedError("Subsampling is not yet implemented")
    max_bins = max_dist // binsize

    # Define chromosomes to scan
    # NOTE: chromosomes smaller than the kernel used for smoothing or the
    # max_dist must be removed.
    # min_size = max((2 * h + 1) * binsize, max_dist)
    # chromlist, chroms_lengths = make_chromlist(mat1, whitelist, blacklist, min_size=min_size)
    # if not len(chromlist):
    #     raise KeyError(
    #         "All chromosomes were too short and have been discarded. " "Try reducing max_dist."
    #     )

    chrom_starts = np.insert(np.cumsum(chroms_lengths), 0, 0)

    # Convert to number of contacts to proportions if needed.
    if subsample is not None:
        if subsample > 1:
            subsample_prop_1 = subsample / mat1.sum()
            subsample_prop_2 = subsample / mat2.sum()
        else:
            subsample_prop_1 = subsample_prop_2 = subsample
        if subsample_prop_1 > 1 or subsample_prop_2 > 1:
            raise ValueError("Subsampling value exceeds matrix contacts")
        if subsample_prop_1 < 0 or subsample_prop_2 < 0:
            raise ValueError("Subsampling values must be positive")
    # Compute SCC values separately for each chromosome
    chroms_scc = np.zeros(len(chroms_lengths))
    for c in range(len(chroms_lengths)):
        # chrom_1 = mat1.matrix(sparse=True, balance=False).fetch(chrom)
        # chrom_2 = mat2.matrix(sparse=True, balance=False).fetch(chrom)
        # slice out chromosome submatrices
        c_start = chrom_starts[c]
        c_end = chrom_starts[c + 1]
        chrom_1 = mat1[c_start:c_end, c_start:c_end]
        chrom_2 = mat2[c_start:c_end, c_start:c_end]
        # Sample 10% contacts and smooth 10 times for this chromosome
        if subsample is not None:
            chrom_1 = mat_process.subsample_contacts(chrom_1, subsample_prop_1)
            chrom_2 = mat_process.subsample_contacts(chrom_2, subsample_prop_2)
        # Trim diagonals which are too far to be scanned to reduce
        # compute time and memory usage
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=sp.SparseEfficiencyWarning)
            chrom_1 = mat_process.diag_trim(chrom_1.todia(), max_bins + h).tocoo()
            chrom_2 = mat_process.diag_trim(chrom_2.todia(), max_bins + h).tocoo()
        smooth_1 = mat_process.smooth(chrom_1, h)
        smooth_2 = mat_process.smooth(chrom_2, h)

        chroms_scc[c] = get_scc(smooth_1, smooth_2, max_bins=max_bins)
    # Compute the genome SCC using the weighted averge of chromosomes
    # SCC by their lengths. NaN values of SCC are not considered
    # This happens when comparing empty diagonals
    nan_scc_mask = ~np.isnan(chroms_scc)
    trunc_scc = chroms_scc[nan_scc_mask]
    trunc_lengths = np.array(chroms_lengths)[nan_scc_mask]
    scc = np.average(trunc_scc, weights=trunc_lengths)
    return scc
