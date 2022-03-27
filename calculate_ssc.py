#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 15:47:59 2022.

@author: mg
"""

from multiprocessing import Pool

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
from tools import scc
from tools.mg_plot import new_fig, set_styling

np.set_printoptions(precision=3)

BIN_SIZE = int(1e5)  # bins are 100 kbps
MAX_DIST = int(5e6)  # max dist to consider is 5 Mbps

H_VAL = 5

N_PARTICLES = 25714  # this is not expected to change anytime soon, so fix it for now


def load_sim_contact_matrix(n_cell: int) -> sp.csr_matrix:
    """Load contact matrix of simulation data, normalised to a maximum of 1.

    Parameters
    ----------
    n_cell : int
        Number of cell for which to load contact matrix

    Returns
    -------
    scipy.sparse.csr_matrix
        CSR sparse contact matrix
    """
    mat = sp.load_npz(f"data/contact_matrices/contact_matrix_cell{n_cell}.npz").tocsr().astype(np.single)

    return mat / mat.max()


# def make_cp_matrix(idxs, N_particles):
#     # idxs, N_particles = args
#     mat = sp.coo_matrix((N_particles, N_particles), dtype="uint8")
#     for i in range(idxs.shape[0]):
#         idx = tuple(idxs[i, :])
#         mat[idx] = 1
#     return mat.tocsr()


def make_hic_contact_matrix(n_cell: int) -> sp.csr_matrix:
    """Make contact matrix from HiC contact pairs.

    Parameters
    ----------
    n_cell : int
        Number of cell for which to make contact matrix

    Returns
    -------
    scipy.sparse.csr_matrix
        CSR sparse contact matrix
    """
    cps = np.unique(
        pd.read_pickle(f"data/contact_pairs/contact_pairs_cell{n_cell}.pkl")[["ind_A", "ind_B"]].to_numpy(), axis=0,
    )

    hic_contact_mat = sp.coo_matrix(
        (np.repeat(1, cps.shape[0]), (cps[:, 0], cps[:, 1])), shape=(N_PARTICLES, N_PARTICLES),
    ).tocsr()

    # cps = cps[cps.argsort(axis=0)[:, 0]]  # sort by first entry

    # make HiC contact matrix
    # hic_contact_mat = np.zeros((N_particles, N_particles), dtype="float32")  # , dtype="int16"

    # for i, j in cps:
    #     hic_contact_mat1[i, j] = 1

    # p = Pool(8)

    # split = np.array_split(cps, 100)

    # arrs = p.starmap(make_cp_matrix, [(s, N_particles) for s in split])

    # hic_contact_mat = np.sum(arrs, axis=0)

    # hic_contact_mat1 /= np.max(hic_contact_mat1)

    # add self-contacts and neighbours
    ex = np.ones(N_PARTICLES)
    dia_data = [ex, ex, ex]
    dia_mat = sp.dia_matrix((dia_data, [-1, 0, 1]), shape=(N_PARTICLES, N_PARTICLES)).tocsr()
    # hic_contact_mat += np.eye(N_particles, k=-1)
    # hic_contact_mat += np.eye(N_particles, k=0)
    # hic_contact_mat += np.eye(N_particles, k=1)

    hic_contact_mat += dia_mat

    assert hic_contact_mat.max() == 1, "HiC contact matrix maximum is not 1"

    return hic_contact_mat


def _f(h):
    print(f"{h}\n", end="")
    mat = np.empty((8, 8))
    for i in range(8):
        for j in range(i + 1):
            n_cell1 = i + 1
            n_cell2 = j + 1

            # hic cross sccs
            # hic_contact_mat1 = make_hic_contact_matrix(n_cell1)
            # hic_contact_mat2 = make_hic_contact_matrix(n_cell2)

            # s = scc.genome_scc(
            #     hic_contact_mat1, hic_contact_mat2, BIN_SIZE, chrom_lengths, MAX_DIST, H_VAL
            # )

            # scc_hic[i, j] = s
            # scc_hic[j, i] = s

            # sim cross sccs
            sim_contact_mat1 = load_sim_contact_matrix(n_cell1)
            sim_contact_mat2 = load_sim_contact_matrix(n_cell2)

            s = scc.genome_scc(sim_contact_mat1, sim_contact_mat2, BIN_SIZE, chroms_lengths, MAX_DIST, h,)

            mat[i, j] = s
            mat[j, i] = s
    return mat


chroms_lengths = pd.read_pickle("data/chrom_lengths.pkl").to_numpy()

# if __name__ == "__main__":

#     # %% h_train

#     print()

#     # n_cell1 = 1
#     # n_cell2 = 2

#     # sim_contact_mat1 = load_sim_contact_matrix(n_cell1)
#     # sim_contact_mat2 = load_sim_contact_matrix(n_cell2)

#     n_cell = 2

#     sim_contact_mat = load_sim_contact_matrix(n_cell)

#     hic_contact_mat = make_hic_contact_matrix(n_cell)

#     h, sccs = scc.h_train(sim_contact_mat, hic_contact_mat, BIN_SIZE, chroms_lengths, MAX_DIST, h_max=10)

#     print(f"Found best h-value {h=}")

#     # h_vals = range(11)

#     # scc_vals = [0.360, 0.133, 0.168, 0.217, 0.259, 0.292, 0.314, 0.330, 0.345, 0.353, 0.355]

#     fig, ax = new_fig()

#     ax.plot(range(len(sccs)), sccs, "C0-o")

#     ax.set_xlabel("h-value")
#     ax.set_ylabel("SCC")

#     ax.set_ylim(0, 0.5)

#     set_styling(ax, x_loc=1)

#     fig.show()

#     # %% Hi-C vs simulated Data

#     print("Starting Hi-C vs Sim SCC")
#     print()

#     sccs = np.empty(8)

#     for n, n_cell in enumerate(range(1, 9)):
#         sim_contact_mat = load_sim_contact_matrix(n_cell)

#         hic_contact_mat = make_hic_contact_matrix(n_cell)

#         s = scc.genome_scc(sim_contact_mat, hic_contact_mat, BIN_SIZE, chroms_lengths, MAX_DIST, h=7)

#         print(f"\tCell {n_cell}: {s:.3}")

#         sccs[n] = s

#     # plt.figure()

#     fig, ax = new_fig()

#     ax.plot(range(1, 9), sccs, "C0o")

#     ax.set_ylim(0, 1)

#     ax.set_xlabel("Cell number")
#     ax.set_ylabel("SCC")

#     set_styling(ax, x_loc=1)

#     fig.show()

# %% Many h-values test

# scc_hic = np.empty((8, 8))

# h_vec = [1, 3, 5, 8, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

# scc_sim = np.empty((len(h_vec), 8, 8))

# for n_h, h in enumerate(h_vec):
#     for i in range(8):
#         for j in range(i + 1):
#             n_cell1 = i + 1
#             n_cell2 = j + 1

#             # hic cross sccs
#             # hic_contact_mat1 = make_hic_contact_matrix(n_cell1)
#             # hic_contact_mat2 = make_hic_contact_matrix(n_cell2)

#             # s = scc.genome_scc(
#             #     hic_contact_mat1, hic_contact_mat2, BIN_SIZE, chrom_lengths, MAX_DIST, H_VAL
#             # )

#             # scc_hic[i, j] = s
#             # scc_hic[j, i] = s

#             # sim cross sccs
#             sim_contact_mat1 = load_sim_contact_matrix(n_cell1)
#             sim_contact_mat2 = load_sim_contact_matrix(n_cell2)

#             s = scc.genome_scc(
#                 sim_contact_mat1,
#                 sim_contact_mat2,
#                 BIN_SIZE,
#                 chroms_lengths,
#                 MAX_DIST,
#                 h,
#             )

#             scc_sim[n_h, i, j] = s
#             scc_sim[n_h, j, i] = s
#         print(
#             f"Progess: {(((i+1) * (i+2)) / 2  + n_h * 8*7) / (len(h_vec)*8*7):.3%}"
#         )
# scc_means = scc_sim.mean(axis=(1, 2))

# p = Pool(8)

# mats = p.map(_f, h_vec)

# scc_means = [mat.mean() for mat in mats]

# np.save("data/scc/scc_hic.npy", scc_hic)
# np.save("data/scc/scc_sim.npy", scc_sim)

# scc_hic = np.load("data/scc/scc_hic.npy")
# scc_sim = np.load("data/scc/scc_sim.npy")

# print()

# print(scc_hic)

# for n_h in range(len(h_vec)):
#     print(scc_sim[n_h, :, :])
#     print()

# plt.plot(h_vec, scc_means, "C0-o")

# plt.figure()

# %%
# print()

print("SCC values between Hi-C data")

print()

scc_mat = np.empty((8, 8))

hic_contact_mats = []

for i in range(8):
    hic_contact_mats += [make_hic_contact_matrix(i + 1)]

for i in range(8):
    for j in range(8):
        scc_mat[i, j] = scc.genome_scc(
            hic_contact_mats[i], hic_contact_mats[j], BIN_SIZE, chroms_lengths, MAX_DIST, h=7
        )

print(scc_mat)

# # %%
# print()

# print("SCC values between simulation data")

# print()

# scc_mat = np.empty((8, 8))

# sim_contact_mats = []

# for i in range(8):
#     sim_contact_mats += [load_sim_contact_matrix(i + 1)]

# for i in range(8):
#     for j in range(8):
#         scc_mat[i, j] = scc.genome_scc(
#             sim_contact_mats[i], sim_contact_mats[j], BIN_SIZE, chroms_lengths, MAX_DIST, h=7
#         )

# print(scc_mat)

# %% Hi-C vs simulated data

sccs = np.empty(8)

for n, n_cell in enumerate(range(1, 9)):
    sim_contact_mat = load_sim_contact_matrix(n_cell)

    hic_contact_mat = make_hic_contact_matrix(n_cell)

    s = scc.genome_scc(sim_contact_mat, hic_contact_mat, BIN_SIZE, chroms_lengths, MAX_DIST, h=7)

    print(f"\tCell {n_cell}: {s:.3}")

    sccs[n] = s

print()

print(f"Min: {sccs.min():.2}")
print(f"Max: {sccs.max():.2}")

# plt.figure()

fig, ax = new_fig()

ax.plot(range(1, 9), sccs, "C0o")

ax.set_ylim(0, 1)

ax.set_xlabel("Cell number")
ax.set_ylabel("SCC")

set_styling(ax, x_loc=1)
