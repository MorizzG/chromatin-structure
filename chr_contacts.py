#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 15 20:33:17 2022.

@author: mg
"""

import gsd.hoomd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse
from tools.mg_plot import new_fig, set_styling

# from multiprocessing import Pool

np.set_printoptions(precision=3)


# %% Load dist and bond data


def calc_contacts(n_cell):
    """
    Calculate the contact ratios between chromosomes for a cell.

    Parameters
    ----------
    n_cell : int
        Cell number

    Returns
    -------
    chr_contact_ratio : numpy.ndarray(shape=(20,20))
        Array of chromosome contact rations

    """
    print(f"Starting cell {n_cell}...")
    print()

    dists = np.load(f"data/dists/dists_cell{n_cell}.npy", allow_pickle=True)

    with gsd.hoomd.open(f"data/trajs/traj_cell{n_cell}.gsd", "rb") as f:

        # pos_all = np.array([f[n].particles.position for n in range(len(f))])

        # pos = pos_all[-1]

        # snap = f[-1]

        all_bonds = f[0].bonds.group
        typeid = f[0].bonds.typeid  # 0 -> bond; 1 -> contact

        # bonds = all_bonds[typeid == 0]
        contacts = all_bonds[typeid == 1]

        # skip the first 5 configurations
        all_pos = np.stack([snap.particles.position for snap in f[5:]])

    # 2x the expected bond length of 1.5
    contact_mat = dists < 3

    return contact_mat


def chr_contact_ratios(contact_mat):
    # chr_lens = pd.read_pickle(f"data/lengths_jan/chr_lens_cell{n_cell}.pkl").to_numpy()
    chr_lens = pd.read_pickle(f"data/chrom_lengths.pkl").to_numpy()

    chr_starts = np.insert(np.cumsum(chr_lens), 0, 0)

    chr_num_contacts = np.array(
        [
            [
                np.count_nonzero(contact_mat[chr_starts[n] : chr_starts[n + 1], chr_starts[m] : chr_starts[m + 1],])
                for n in range(20)
            ]
            for m in range(20)
        ]
    )

    chr_contact_ratio = chr_num_contacts / np.outer(chr_lens, chr_lens)

    # set self-contact ratio to zero since it distorts the scale too much
    # other contact-ratios are basically zero compared to self-contacts
    for i in range(20):
        # for j in range(i, 20):
        #     chr_contact_ratio[i, j] = 0
        chr_contact_ratio[i, i] = 0

    return chr_contact_ratio

    # chr_contact_ratios[n_cell] = chr_contact_ratio


# contact_mat = calc_contacts(1)

# %% Calculate number of contacts captured by Hi-C

# for n_cell in range(1, 9):
#     print(f"Cell {n_cell}")

#     print()

#     contact_mat_sparse = scipy.sparse.load_npz(f"data/contact_matrices/contact_matrix_cell{n_cell}.npz")

#     contact_mat = contact_mat_sparse.toarray()

#     num_neighbour_contacts = 0

#     for i in np.arange(contact_mat.shape[0] - 1):
#         num_neighbour_contacts += contact_mat[i, i]
#         num_neighbour_contacts += contact_mat[i, i + 1]
#         num_neighbour_contacts += contact_mat[i + 1, i]

#     num_neighbour_contacts += contact_mat[-1, -1]

#     # count number of contacts, remove diagonal (self-contacts) and divide by two
#     # since each contact is twice in the symmetric matrix
#     # divide by 100 since sum over 100 frames
#     N_all_contacts = (np.sum(contact_mat) - num_neighbour_contacts) // 2 / 100

#     with gsd.hoomd.open(f"data/trajs/traj_cell{n_cell}.gsd", "rb") as f:

#         # pos_all = np.array([f[n].particles.position for n in range(len(f))])

#         # pos = pos_all[-1]

#         # snap = f[-1]

#         all_bonds = f[0].bonds.group
#         typeid = f[0].bonds.typeid  # 0 -> bond; 1 -> contact

#         # bonds = all_bonds[typeid == 0]
#         contacts = all_bonds[typeid == 1]

#     N_contacts = len(contacts)

#     print(f"\tContacts prespecified: {N_contacts}")

#     print(f"\tContacts simulation: {N_all_contacts:.0f}")

#     print(f"\tPercentage of contacts specified: {N_contacts / N_all_contacts:.1%}")

#     print()


# %% Contact strength between chromosomes

n_cell = 1

fig, axes = new_fig(nrows=4, ncols=2)  # , sharex=True, sharey=True

axes = axes.flatten()

for n_cell in range(1, 9):

    contact_mat_sparse = scipy.sparse.load_npz(f"data/contact_matrices/contact_matrix_cell{n_cell}.npz")

    contact_mat = contact_mat_sparse.toarray()

    chr_contact_mat = chr_contact_ratios(contact_mat)

    # print(chr_contact_ratio)

    ax = axes[n_cell - 1]

    ax.imshow(chr_contact_mat, cmap="Purples")

    ax.set_xticks([])
    ax.set_yticks([])

plt.subplots_adjust(left=0.1, right=0.8, bottom=0.05, top=0.95, wspace=0.0, hspace=0.0)

axes[6].set_xlabel("Chr 0..19, X")
axes[0].set_ylabel("Chr 0..19, X")

cax = plt.axes([0.82, 0.05, 0.05, 0.9])

fig.colorbar(mpl.cm.ScalarMappable(norm=None, cmap="Purples"), cax=cax)

# ax.set_xticks(ticks=None)  #, labels=None)

# %% OLD

# X = np.arange(contact_mat.shape[0], step=100)

# chrom_contact_mat = np.array(
#     [[np.sum(contact_mat[X[i] : X[i + 1], X[j] : X[j + 1]]) for i in range(len(X) - 1)] for j in range(len(X) - 1)]
# )

# fig, ax = plt.subplots(figsize=(8, 6))

# im = ax.imshow(chrom_contact_mat)
# plt.colorbar(im)

# p = Pool()

# # chr_contact_ratios = list(map(calc_contact_ratios, range(1, 9)))

# # %% Display chromosome contact ratios graphically

# print("Plotting...")

# chr_contact_ratios /= np.max(chr_contact_ratios)

# norm = mpl.colors.Normalize(vmin=0, vmax=np.max(chr_contact_ratios))

# fig, axes = plt.subplots(3, 3, figsize=(12, 14))

# # axes = axes.flatten()

# axes_flat = axes.ravel()

# axes_flat[-1].set_frame_on(False)
# axes_flat[-1].tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

# for n_cell in range(8):
#     axes_flat[n_cell].imshow(chr_contact_ratios[n_cell], norm=norm, cmap="plasma")

#     # axes_flat[n_cell].tick_params(
#     #     left=False, bottom=False  # , labelleft=False, labelbottom=False
#     # )

# fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="plasma"), ax=axes)

# # set_styling(ax)


# # %% Calc means and stds

# ratios = np.stack(chr_contact_ratios)

# ratios = np.moveaxis(ratios, 0, -1)

# means = ratios.mean(axis=-1)

# stds = ratios.std(axis=-1)

# fig, axes = plt.subplots(1, 2, figsize=(14, 6))


# im1 = axes[0].imshow(means, cmap="plasma")
# fig.colorbar(im1, cmap="plasma", ax=axes[0])

# im2 = axes[1].imshow(stds, cmap="plasma")
# fig.colorbar(im2, cmap="plasma", ax=axes[1])
