#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 15 19:45:29 2022.

@author: mg
"""


# %% Imports

import argparse

import gsd.hoomd
import numpy as np

# import pandas as pd
import seaborn as sns
from scipy.stats import rv_histogram
from tools.mg_plot import new_fig, set_styling

# import matplotlib.pyplot as plt


# from multiprocessing import Pool


# from scipy.signal import savgol_filter


# %% Parse args


# parser = argparse.ArgumentParser(description="Calculate distance distributions for bonds and contacts")

# arg_group = parser.add_argument(
#     "cell_n", action="store", nargs="?", type=int, default=1, help="Which cell to calculate",  # , nargs="*"
# )

# args = parser.parse_args()

# cell_n = args.cell_n

n_cell = 2


# %% Read trajectory data

# for n_cell in range(1, 9):

print(f"Cell {n_cell}:")

print()

with gsd.hoomd.open(f"data/trajs/traj_cell{n_cell}.gsd", "rb") as f:

    # pos_all = np.array([f[n].particles.position for n in range(len(f))])

    # pos = pos_all[-1]

    # snap = f[-1]

    all_bonds = f[0].bonds.group
    typeid = f[0].bonds.typeid  # 0 -> bond; 1 -> contact

    bonds = all_bonds[typeid == 0]
    contacts = all_bonds[typeid == 1]

    # skip the first 5 configurations
    all_pos = np.stack([snap.particles.position for snap in f[5:]])  # 5:

# %% Read distance data

# all_dists = np.load(f"data/dists/dists_cell{cell_n}.npy", allow_pickle=True)

# all_dists = np.fromiter(
#     (all_dists[n, m] for n in range(len(all_dists)) for m in range(n)),
#     dtype="float32",
#     count=(len(all_dists) * (len(all_dists) - 1) // 2),
# )

# %% Calculate bond and contact distance distribution

bond_dists = np.linalg.norm(np.concatenate([pos[bonds[:, 0]] - pos[bonds[:, 1]] for pos in all_pos]), axis=1)

contact_dists = np.linalg.norm(
    np.concatenate([pos[contacts[:, 0]] - pos[contacts[:, 1]] for pos in all_pos]), axis=1
)

rv_bonds = rv_histogram(np.histogram(bond_dists, "auto", density=True))

rv_contacts = rv_histogram(np.histogram(contact_dists, "auto", density=True))

# rv_all = rv_histogram(np.histogram(np.ravel(all_dists), "auto", density=True))

# %% Print averages and quartiles

print(f"\tBonds average distance: {rv_bonds.mean():.2f}")
print(f"\tBonds 99.73% quartile: {rv_bonds.ppf(.9973):.2f}")

print()

print(f"\tContacts average distance: {rv_contacts.mean():.2f}")
print(f"\tContacts 99.73% quartile: {rv_contacts.ppf(.9973):.2f}")

print()


# %% Plot contact and bond distance distribution

# fig, ax = new_fig()

# X = np.linspace(0, 3, num=1000)


# def pdf_bonds_smooth(x, window_length=21, window_width=0.07):
#     X = np.linspace(x - window_width, x + window_width, num=window_length)
#     return savgol_filter(rv_bonds.pdf(X), window_length, 3, axis=0)[10, :]


# def pdf_contacts_smooth(x, window_length=21, window_width=0.1):
#     X = np.linspace(x - window_width, x + window_width, num=window_length)
#     return savgol_filter(rv_contacts.pdf(X), window_length, 3, axis=0)[10, :]


# ax.plot(X, rv_bonds.pdf(X), "C0-", label="bond length PDF")
# ax.plot(X, rv_contacts.pdf(X), "C1-", label="contact length PDF")

# (ylim_low, ylim_high) = ax.get_ylim()

# ax.vlines(bond_dists.mean(), ylim_low - 1, ylim_high + 1, colors="C0")
# ax.vlines(contact_dists.mean(), ylim_low - 1, ylim_high + 1, colors="C1")

# ax.set_xlim(0, 3)
# ax.set_ylim(ylim_low, ylim_high)

# ax.set_xlabel("distance")
# ax.set_ylabel("PDF")

# set_styling(ax)

# ax.legend()


# ax2 = ax.twinx()

# # X = np.linspace(0, 5, num=1000)

# (line3,) = ax2.plot(X, rv_all.pdf(X), "C2-", label="all dists PDF")

# lines = [line1, line2, line3]

# ax.legend(lines, [line.get_label() for line in lines])


# %% Plot contact and bond distribution using seaborn

fig, ax = new_fig()

sns.kdeplot(bond_dists, ax=ax, label="bond length PDF")
sns.kdeplot(contact_dists, ax=ax, label="contact length PDF")

# (ylim_low, ylim_high) = ax.get_ylim()

# ax.vlines(bond_dists.mean(), ylim_low - 1, ylim_high + 1, colors="C0")
# ax.vlines(contact_dists.mean(), ylim_low - 1, ylim_high + 1, colors="C1")

ax.axvline(bond_dists.mean(), color="C0")
ax.axvline(contact_dists.mean(), color="C1")

ax.set_xlim(0.5, 3)
# ax.set_ylim(ylim_low, ylim_high)

ax.set_xlabel("distance")
ax.set_ylabel("PDF")

set_styling(ax)

ax.legend(loc=(0.6, 0.83))
