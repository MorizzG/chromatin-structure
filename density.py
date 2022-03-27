# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 21:30:54 2022.

@author: mg
"""
import gsd.hoomd
import numpy as np


def min_max_cube(a):
    return tuple((np.min(a[:, i]), np.max(a[:, i])) for i in [0, 1, 2])


def density_grid(pos):
    lims = min_max_cube(pos)

    lims = [(np.floor(lower), np.ceil(higher)) for (lower, higher) in lims]

    side_lens = [int(higher - lower) for (lower, higher) in lims]


n_cell = 2

traj = gsd.hoomd.open(f"data/trajs/traj_cell{n_cell}.gsd")

pos_all = np.stack([snap.particles.position for snap in traj], axis=0)
