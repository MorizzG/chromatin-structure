#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 29 18:30:07 2022

@author: mg
"""

# %% Imports and globals

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.signal
import scipy.sparse

n_cell = 1

# %% Load data

df_contact_pairs = pd.read_pickle(f"data/contact_pairs_jan/contact_pairs_cell{n_cell}.pkl")

lengths = pd.read_pickle(f"data/lengths_jan/chr_lens_cell{n_cell}.pkl")


# %% Make contact matrix

cp = df_contact_pairs[["ind_A", "ind_B"]].values

N = np.sum(lengths)

mat = np.zeros((N, N), dtype="int8")

for (i, j) in cp:
    mat[i, j] = 1
    mat[j, i] = 1

# w = scipy.signal.get_window("hamming", 11)

# w = w / np.max(w)

# w = np.outer(w, w)

# mat = scipy.signal.convolve(mat, w, mode="same")

X = np.arange(N, step=100)

mat_thinned = np.array(
    [[np.sum(mat[X[i] : X[i + 1], X[j] : X[j + 1]]) for i in range(len(X) - 1)] for j in range(len(X) - 1)]
)


# %% Show contact matrix

fig, ax = plt.subplots(figsize=(8, 6))

im = ax.imshow(mat_thinned)
plt.colorbar(im)
