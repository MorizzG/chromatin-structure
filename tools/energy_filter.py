# -*- coding: utf-8 -*-
"""
Created on Tue Feb 22 19:28:54 2022.

@author: mg
"""

import numpy as np

# import gsd.hoomd
import pandas as pd

# CUTOFF_PERCENTAGE = 1


def filter_energy(traj, n_cell, n_chrom=None, cutoff_percentage=1):
    log_path = f"data/logs/log_cell{n_cell}" + (f"_chrom{n_chrom}" if n_chrom else "") + ".log"
    df_log = pd.read_csv(log_path, sep="\t")

    energies = df_log[["potential_energy"]].to_numpy()

    min_energy = energies.min()

    low_energy_idxs = np.nonzero(energies < (100 + cutoff_percentage) / 100 * min_energy)[0]

    return [traj[int(idx)] for idx in low_energy_idxs]
