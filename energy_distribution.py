#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 15 19:45:29 2022.

@author: mg
"""

# import argparse

import numpy as np
import pandas as pd
from matplotlib.ticker import FormatStrFormatter
from tools.mg_plot import new_fig, set_styling

# import matplotlib.pyplot as plt


# %% Parse args


# parser = argparse.ArgumentParser(description="Display potential energy distribution for simulation")

# arg_group = parser.add_argument(
#     "n_cell", action="store", nargs="?", type=int, default=1, help="Which cell to calculate",  # , nargs="*"
# )

# args = parser.parse_args()

# n_cell = args.n_cell

n_cell = 8

n_chrom = 1

ext = ""
# ext = f"_chrom{n_chrom}"


# %% Plot energy distribution

# for n_cell in range(1, 9):

df_log = pd.read_csv(f"data/logs/log_cell{n_cell}{ext}.log", sep="\t")  # _chrom5

energies = df_log["potential_energy"]

print(f"Cell {n_cell}")

print()

print(f"\tAverage percentage deviation: {energies[5:].std() / energies[5:].mean():.2%}")

print()

fig, ax = new_fig()

ax.plot(energies, "C0.")

# idxs = np.array([ 2,  6, 22, 49, 64, 74, 85, 90, 96])
# ax.plot(idxs, energies[idxs], "C3o", markersize=3)

ax.set_xlabel("frame")
ax.set_ylabel("potential energy")

# ax.xaxis.set_major_formatter(FormatStrFormatter('% 1.2f'))
# ax.ticklabel_format(axis="y", style="plain")

set_styling(ax)  # , y_loc=(2e5, 1e5)


print(f"Lowest energy frame: {energies.argmin()}")
