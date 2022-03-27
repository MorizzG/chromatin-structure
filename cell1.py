# -*- coding: utf-8 -*-
"""
Created on Mon Feb 21 14:29:41 2022.

@author: mg
"""
# %% Imports
import gsd.hoomd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from tools import rmsd
from tools.mg_plot import new_fig, set_styling

# %% Globs

ENERGY_CUTOFF = 1.495e7

# %% Load data

n_cell = 1

ext = ""

df_log = pd.read_csv(f"data/logs/log_cell{n_cell}{ext}.log", sep="\t")
f = gsd.hoomd.open(f"data/trajs/traj_cell{n_cell}{ext}.gsd")

energies = df_log["potential_energy"].to_numpy()

low_energy_configs = np.nonzero(energies < ENERGY_CUTOFF)[0]

# %% Plot energies

fig, ax = new_fig()

ax.plot(energies, "C0.")
ax.axhline(ENERGY_CUTOFF, color="C1", label="Cutoff energy")

ax.set_xlabel("frame")
ax.set_ylabel("potential energy")

ax.legend(loc=(0.65, 0.85))
set_styling(ax)

# %% RMSD

pos_all = np.stack([snap.particles.position for snap in f], axis=0)

rmsds = np.empty(105)

ref_idx = low_energy_configs[-1]

for i, pos in enumerate(pos_all):
    rmsds[i] = rmsd.rmsd(pos, pos_all[ref_idx, :, :])

rmsds[ref_idx] = 0

# %% Plot RMSDs

fig, ax = new_fig()

ax.plot(rmsds, "C0-")

ax.plot(low_energy_configs, rmsds[low_energy_configs], "C3o", markersize=3, label="Low energy configs")

ax.set_xlabel("frame")
ax.set_ylabel("RMSD")

ax.set_ylim(top=11.0)

ax.legend(loc=(0.60, 0.88))
set_styling(ax)
