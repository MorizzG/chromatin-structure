#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 19:15:36 2022

@author: mg
"""

# %% Imports

import numpy as np
import pandas as pd
import scipy.sparse

# %% Global functions and constants

n_cell1 = 1
n_cell2 = 2


def mat_norm(A):
    return np.sum(np.abs(A))


def unit(A):
    return A / mat_norm(A)


def d(A, B):
    return mat_norm(A - B) / (mat_norm(A) + mat_norm(B))
    # return mat_norm(unit(A) - unit(B)) / 2


# %% Load data

# since the contact matrices are sums of 100 configurations, normalise by dividing by 100
cell1_sim = (
    scipy.sparse.load_npz(f"data/contact_matrices/contact_matrix_cell{n_cell1}.npz").toarray().astype("float32")
) / 100
cell2_sim = (
    scipy.sparse.load_npz(f"data/contact_matrices/contact_matrix_cell{n_cell2}.npz").toarray().astype("float32")
) / 100

assert cell1_sim.shape[0] == cell1_sim.shape[1] == cell2_sim.shape[0] == cell2_sim.shape[1]
N_particles = cell1_sim.shape[0]

cps1 = pd.read_pickle(f"data/contact_pairs/contact_pairs_cell{n_cell1}.pkl")[["ind_A", "ind_B"]].to_numpy()
cps2 = pd.read_pickle(f"data/contact_pairs/contact_pairs_cell{n_cell2}.pkl")[["ind_A", "ind_B"]].to_numpy()

# %% Make cp data cells


cell1_cp = np.zeros((N_particles, N_particles), dtype="float32")  # , dtype="int16"
cell2_cp = np.zeros((N_particles, N_particles), dtype="float32")  # , dtype="int16"

for i, j in cps1:
    cell1_cp[i, j] += 1

for i, j in cps2:
    cell2_cp[i, j] += 1

cell1_cp /= np.max(cell1_cp)
cell2_cp /= np.max(cell2_cp)

# we count chromosome boundaries as contact as well
# but this only introduces a neglegible error
# for i in range(N_particles - 1):
#     cell1_cp[i, i] = 1
#     cell1_cp[i, i + 1] = 1
#     cell1_cp[i + 1, i] = 1

#     cell2_cp[i, i] = 1
#     cell2_cp[i, i + 1] = 1
#     cell2_cp[i + 1, i] = 1

cell1_cp += np.eye(N_particles, k=-1)
cell1_cp += np.eye(N_particles, k=0)
cell1_cp += np.eye(N_particles, k=1)

cell2_cp += np.eye(N_particles, k=-1)
cell2_cp += np.eye(N_particles, k=0)
cell2_cp += np.eye(N_particles, k=1)

cell1_cp[-1, -1] = 1
cell1_cp[-1, -1] = 1

# %% Calculate

a = d(cell1_sim, cell2_sim)

print(a)
