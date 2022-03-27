# -*- coding: utf-8 -*-
"""
Created on Thu Mar 10 18:59:25 2022.

@author: mg
"""

from multiprocessing import Pool
from os import path

import gsd.hoomd
import numpy as np
from tools.rmsd import rmsd


def _rmsd_worker(n_cell, ref_idx):
    return cell_file_rmsds(f"../data/trajs/traj_cell{n_cell}.gsd", ref_idx=ref_idx)


def create_all_rmsd_files():
    """Calculate cell RMSDs with respect to each frame and save to file for faster analysis."""
    for n_cell in range(1, 9):
        try:
            rmsds = np.load(f"../data/rmsds/rmsds_cell{n_cell}.npy")
        except FileNotFoundError:
            print(f"RMSD file not found for cell {n_cell}, rerunning calculation")

            pool = Pool()

            arr = pool.starmap(_rmsd_worker, zip(105 * [n_cell], range(105)))

            rmsds = np.stack(arr)

            np.save(f"../data/rmsds/rmsds_cell{n_cell}.npy", rmsds)


# %% Parse args


# %% Load traj and calculate rmsd


def cell_file_rmsds(filepath, ref_idx=-1):
    """
    Calculate RMSD of all configurations in a cell with respect to a reference configuration.

    Parameters
    ----------
    filepath : str
        file path to gsd file of cell
    ref_frame : int, optional
        Reference frame The default is -1.

    Returns
    -------
    numpy.array
        Array of RMSDs

    """
    if not path.exists(filepath):
        raise ValueError(f"Error: file {filepath} does not exist")
    # all_pos = np.load(f"data/trajs_np/traj_np_cell{cell_n}.npy")
    traj = gsd.hoomd.open(filepath)

    all_pos = np.stack([snap.particles.position for snap in traj], axis=0)

    pos_ref = all_pos[ref_idx, :, :]

    rmsds = np.empty(all_pos.shape[0])

    for idx in range(all_pos.shape[0]):
        pos = all_pos[idx, :, :]

        # mean_pos = np.mean(pos, axis=0)

        # norm_pos = np.linalg.norm(pos - mean_pos)

        # mtx1, mtx2, M = scipy.spatial.procrustes(pos, pos_ref)

        # rmsd[idx] = norm_pos * np.sqrt(3 * ((mtx1 - mtx2) ** 2).mean())
        rmsds[idx] = rmsd(pos, pos_ref)
    rmsds[ref_idx] = 0

    return rmsds


if __name__ == "__main__":
    create_all_rmsd_files()
