#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 15:56:38 2022.

@author: mg
"""

import argparse

# import sys
from multiprocessing import Pool

# import pandas as pd
import gsd.hoomd
import numpy as np

# import scipy.spatial


NUM_THREADS = 8


# %% Parse args

parser = argparse.ArgumentParser(description="Calculate chromosome contactivity")

arg_group = parser.add_argument(
    "n_cell", action="store", nargs="?", type=int, default=1, help="Which cell to calculate",
)

args = parser.parse_args()

n_cell = args.n_cell


def calc_all_dists(n):
    return np.linalg.norm(np.array([pos[n] - pos[m] for m in range(len_pos)]), axis=1)


f = gsd.hoomd.open(f"data/trajs/traj_cell{n_cell}.gsd", "rb")

snap = f[-1]

pos = snap.particles.position

len_pos = len(pos)


if __name__ == "__main__":
    print("starting")

    p = Pool(NUM_THREADS)

    dists = np.stack(p.map(calc_all_dists, range(len_pos)))

    # dists = scipy.spatial.distance_matrix(pos, pos)

    np.save(f"data/dists/dists_cell{n_cell}", dists)
