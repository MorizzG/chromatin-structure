# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 16:09:55 2022.

@author: mg
"""
import warnings

# %% Import
import gsd.hoomd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as ss
from tools import mat_process, rmsd, scc
from tools.mg_plot import new_fig, set_styling

np.set_printoptions(precision=3)

# %% Loading data

n_cell = 8
n_chrom = 1

mean_rmsds_to_cell = np.empty(8)
int_contacts = np.empty(8)

for n_cell in range(1, 9):

    f_cell = gsd.hoomd.open(f"data/trajs/traj_cell{n_cell}.gsd")
    f_chrom = gsd.hoomd.open(f"data/trajs/traj_cell{n_cell}_chrom{n_chrom}.gsd")

    chroms_lens = pd.read_pickle("data/chrom_lengths.pkl")
    chroms_lens_cum = np.cumsum(chroms_lens)

    # %% Data wrangling

    pos_all_cell = np.stack([snap.particles.position for snap in f_cell])

    pos_all_chrom = np.stack([snap.particles.position for snap in f_chrom])

    chrom_start = chroms_lens_cum[n_chrom - 1] if n_chrom > 1 else 0
    chrom_end = chroms_lens_cum[n_chrom]

    # filter out positions of selected chromosomel
    pos_all_cell = pos_all_cell[:, chrom_start:chrom_end, :]

    assert (
        pos_all_cell.shape == pos_all_chrom.shape
    ), f"pos_all_cell shape {pos_all_cell.shape} and pos_all_chrom shape {pos_all_chrom.shape} don't match"

    # %% Cell average trajectory

    # pos_all = np.stack([snap.particles.position for snap in f_cell])

    pos_filtered = None

    idx_filtered = None

    if n_cell == 1:
        # filter frames by energy, drop high energy frames
        ENERGY_CUTOFF = 1.495e7  # from cell1.py

        df_log = pd.read_csv("data/logs/log_cell1.log", sep="\t")
        energies = df_log["potential_energy"].to_numpy()

        low_energy_configs = np.nonzero(energies < ENERGY_CUTOFF)[0]

        pos_filtered = pos_all_cell[low_energy_configs, :, :]
    elif n_cell == 5:
        # drop all frames before 39, since high energy
        pos_filtered = pos_all_cell[39:, :, :]
    # elif n_cell == 4:
    #     drop_idxs = np.array([2, 6, 22, 49, 64, 74, 85, 90, 96])
    #     pos_filtered = np.delete(pos_all, drop_idxs, axis=0)
    else:
        # for all other cells, simply drop first 5 frames
        pos_filtered = pos_all_cell[5:, :, :]

    cell_avg_traj = rmsd.average_trajectory(pos_filtered)

    # %% RMSDS

    rmsds_cell = np.empty(105)

    for i in range(104):
        rmsds_cell[i] = rmsd.rmsd(pos_all_cell[i, :, :], cell_avg_traj)

    rmsds_chrom = np.empty(105)
    df_log = pd.read_csv(f"data/logs/log_cell{n_cell}_chrom{n_chrom}.log", sep="\t")  # _chrom5
    energies = df_log["potential_energy"]
    chrom_ref_idx = np.argmin(energies)

    for i in range(105):
        if i == chrom_ref_idx:
            continue

        rmsds_chrom[i] = rmsd.rmsd(pos_all_chrom[i], pos_all_chrom[chrom_ref_idx])

    rmsds_chrom[chrom_ref_idx] = 0

    cross_rmsds = np.empty(105)

    for i_chrom in range(105):
        r = rmsd.rmsd(cell_avg_traj, pos_all_chrom[i_chrom])
        cross_rmsds[i_chrom] = r

    # %% Plot RMSDs

    # fig, ax = new_fig()

    # ax.plot(rmsds_cell)

    # ax.set_xlabel("Configuration")
    # ax.set_ylabel("RMSD")
    # set_styling(ax)

    # fig.suptitle("Cell RMSD")

    # fig, ax = new_fig()

    # ax.plot(rmsds_chrom)

    # ax.set_xlabel("Configuration")
    # ax.set_ylabel("RMSD")
    # set_styling(ax)

    # fig.suptitle("Chrom RMSD")

    # fig, ax = new_fig()

    # im = ax.plot(cross_rmsds)

    # ax.set_xlabel("Chrom Configuration")
    # ax.set_ylabel("Cell Configuration")
    # set_styling(ax)

    # fig.suptitle("Cell vs Chrom")

    print(f"Cell {n_cell}")

    # print(f"\tAverage chrom RMSD:         {rmsds_chrom.mean():.1f} +/- {rmsds_chrom.std():.1f}")

    # print(f"\tAverage cell vs chrom RMSD: {cross_rmsds.mean():.1f} +/- {cross_rmsds.std():.1f}")

    mean_rmsds_to_cell[n_cell - 1] = rmsds_chrom.mean()

    # print()

    # fig, ax = new_fig()

    # ax.plot(cross_rmsds[104, :])

    # ax.set_xlabel("Configuration")
    # ax.set_ylabel("RMSD")
    # set_styling(ax)

    # fig.suptitle("1 Cell vs all Chrom")

    # %% Internal vs external contact ratio in simulation

    # contact_mat = sp.load_npz(f"data/contact_matrices/contact_matrix_cell{n_cell}.npz").toarray()

    # num_int_contacts = contact_mat[chrom_start:chrom_end, chrom_start:chrom_end].sum()

    # num_ext_contacts = (
    #     contact_mat[chrom_start:chrom_end, :chrom_start].sum() + contact_mat[chrom_start:chrom_end, chrom_end:].sum()
    # )

    # print(f"\tInternal contacts: {num_int_contacts}")
    # print(f"\tExternal contacts: {num_ext_contacts}")
    # print(f"\tPercentage internal contacts: {num_int_contacts / (num_int_contacts + num_ext_contacts):.1%}")

    # print()

    # %% Internal vs external contacts in Hi-C

    df_contact_pairs = pd.read_pickle(f"data/contact_pairs/contact_pairs_cell{n_cell}.pkl")

    df_int = df_contact_pairs[(df_contact_pairs["chr_A"] == n_chrom) & (df_contact_pairs["chr_B"] == n_chrom)]
    num_int_contacts = np.unique(df_int.to_numpy(), axis=0).shape[0]

    df_ext = df_contact_pairs[(df_contact_pairs["chr_A"] == n_chrom) & (df_contact_pairs["chr_B"] != n_chrom)]
    num_ext_contacts = np.unique(df_ext.to_numpy(), axis=0).shape[0]

    df_ext = df_contact_pairs[(df_contact_pairs["chr_A"] != n_chrom) & (df_contact_pairs["chr_B"] == n_chrom)]
    num_ext_contacts += np.unique(df_ext.to_numpy(), axis=0).shape[0]

    print(f"\tInternal contacts: {num_int_contacts}")
    print(f"\tExternal contacts: {num_ext_contacts}")
    print(f"\tPercentage internal contacts: {num_int_contacts / (num_int_contacts + num_ext_contacts):.0%}")

    print()

    int_contacts[n_cell - 1] = num_int_contacts

fig, ax = new_fig()

ax.plot(int_contacts, mean_rmsds_to_cell, "C0o")

for n_cell, (x, y) in enumerate(zip(int_contacts, mean_rmsds_to_cell)):
    ax.text(x, y, f"Cell {n_cell+1}")

ax.set_xlabel("Number of internal contacts")
ax.set_ylabel("Mean RMSD of invididual chromosome")

print(f"Pearson R^2: {ss.pearsonr(int_contacts, mean_rmsds_to_cell)[0]:.2f}")


# %% RMSD between chroms


# %% SCC between chroms


def load_sim_contact_matrix(n_cell, n_chrom) -> sp.csr_matrix:
    mat = sp.load_npz(f"data/contact_matrices/contact_matrix_cell{n_cell}_chrom{n_chrom}.npz").tocsr().astype(np.single)

    return mat / mat.max()


print("SCC values between simulation data")

print()

scc_mat = np.empty((8, 8))

h = 7
max_bins = int(5e6) // int(1e5)

sim_contact_mats = []

for n_cell in range(1, 9):
    chrom = load_sim_contact_matrix(n_cell, n_chrom)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=sp.SparseEfficiencyWarning)
        chrom = mat_process.diag_trim(chrom.todia(), max_bins + h).tocoo()
    chrom = mat_process.smooth(chrom, h)
    sim_contact_mats += [chrom]

for i in range(8):
    for j in range(8):
        scc_mat[i, j] = scc.get_scc(sim_contact_mats[i], sim_contact_mats[j], max_bins)

print(scc_mat)
