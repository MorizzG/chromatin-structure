#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 29 21:09:32 2022

@author: mg
"""

# %% Imports

import argparse
import sys
from multiprocessing import Pool

import gsd.hoomd
import numpy as np
import scipy.sparse
import scipy.spatial

# import pandas as pd


# import matplotlib as mpl
# import matplotlib.pyplot as plt


# from tools.mg_plot import new_fig, set_styling


NUM_THREADS = 4


# %% Parse args

parser = argparse.ArgumentParser(description="Calculate chromosome contactivity")

parser.add_argument(
    "n_cell", action="store", type=int, help="Which cell to calculate",
)

parser.add_argument(
    "n_chrom", action="store", nargs="?", type=int, default=None, help="Which cell to calculate",
)

args = parser.parse_args()

n_cell = args.n_cell

n_chrom = args.n_chrom


# %% Read trajectory data

ext = f"_chrom{n_chrom}" if n_chrom else ""

f = gsd.hoomd.open(f"data/trajs/traj_cell{n_cell}{ext}.gsd", "rb")

pos_all = np.stack([snap.particles.position for snap in f[5:]], axis=0)

N_particles = pos_all.shape[1]


# %% Calculate contact matrix


# def calc_contact_matrix(n):
#     print(f"Calculating snap {n}\n", end="")
#     pos = pos_all[n, :, :]
#     dists = scipy.spatial.distance_matrix(pos, pos)

#     return dists < 3


# num_contacts = np.zeros((N_particles, N_particles), dtype="int8")

# snap_idxs = range(pos_all.shape[0])

# p = Pool(NUM_THREADS)

# while snap_idxs:
#     idxs = snap_idxs[:NUM_THREADS]
#     snap_idxs = snap_idxs[NUM_THREADS:]

#     contacts_list = p.map(calc_contact_matrix, idxs)

#     # for contacts in contacts_list:
#     #     num_contacts += contacts

#     num_contacts += np.sum(contacts_list, axis=0)


# def count_contacts(args):
#     (i, j) = args
#     # print(f"Calculating pos {i,j}\n", end="")

#     dists = np.linalg.norm(pos_all[:, i, :] - pos_all[:, j, :], axis=1)
#     return np.sum(dists < 3)


def count_contacts(i):
    row_rep = np.repeat(pos_all[:, i, :][:, None, :], N_particles, axis=1)
    dists = np.linalg.norm(row_rep - pos_all, axis=-1)
    return np.sum(dists < 3, axis=0)


if __name__ == "__main__":

    p = Pool(NUM_THREADS)

    # num_contacts = np.array(
    #     p.map(
    #         count_contacts, ((i, j) for i in range(N_particles) for j in range(N_particles))
    #     )
    # ).reshape((N_particles, N_particles))

    num_contacts = np.empty((N_particles, N_particles))

    # for i in range(N_particles):
    #     print(f"Doing row {i}")
    #     for j in range(N_particles):
    #         num_contacts[i, j] = count_contacts((i, j))

    # for i in range(N_particles):
    #     print(f"Doing row {i}")
    #     num_contacts[i, :] = count_contacts(i)

    I_IDX = range(N_particles)

    print(f"{0:.1%} done")

    while I_IDX:
        idxs = I_IDX[: 100 * NUM_THREADS]
        I_IDX = I_IDX[100 * NUM_THREADS :]

        ccs = p.map(count_contacts, idxs)

        for (i, cc) in zip(idxs, ccs):
            num_contacts[i] = cc

        print(f"{idxs[-1] / N_particles:.1%} done")

    num_contacts_sparse = scipy.sparse.bsr_matrix(num_contacts, dtype="uint8")

    scipy.sparse.save_npz(f"data/contact_matrices/contact_matrix_cell{n_cell}{ext}", num_contacts_sparse)
