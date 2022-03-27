#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 20:12:34 2022.

@author: mg
"""
import gsd.hoomd
import matplotlib.pyplot as plt
import numpy as np
import tools.rmsd as rmsd
from tools import scc
from tools.constants import NUM_PARTICLES
from tools.mg_plot import new_fig, set_styling

# %% Setup

np.set_printoptions(precision=3)


# %% Parse args


# parser = argparse.ArgumentParser(description="Calculate chromosome contactivity")

# arg_group = parser.add_argument(
#     "cell_n",
#     action="store",
#     nargs="?",
#     type=int,
#     default=1,
#     help="Which cell to calculate",
# )

# args = parser.parse_args()

# cell_n = 2  # args.cell_n


# %% Plot RMSD for trajectory

n_cell = 3

n_chrom = 1

ext = f"_chrom{n_chrom}"

# ext = ""

# for n_cell in range(1, 9):

traj = gsd.hoomd.open(f"data/trajs/traj_cell{n_cell}{ext}.gsd")

pos_all = np.stack([snap.particles.position for snap in traj])

# pos_avg = rmsd.average_trajectory(pos_all)

rmsds = np.empty(pos_all.shape[0])

ref_idx = 22

for n_frame in range(pos_all.shape[0]):
    # rmsds[n_frame] = rmsd.rmsd(pos_all[n_frame, :, :], pos_avg)
    rmsds[n_frame] = rmsd.rmsd(pos_all[n_frame,:,:], pos_all[ref_idx,:,:])

rmsds[ref_idx] = 0

fig, ax = new_fig()

ax.plot(rmsds)

ax.set_xlabel("frame")
ax.set_ylabel("RMSD")

set_styling(ax)

print()

print(f"Cell {n_cell} mean RMSD: {rmsds.mean():.1f} +/- {rmsds.std():.1f}")


# %% Plot RMSD to last config for all cells

# fig, axes = plt.subplots(3, 3, figsize=(10, 10), sharey=True)

# axes_flat = axes.ravel()

# axes_flat[-1].set_frame_on(False)
# axes_flat[-1].tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)

# for i, ax in enumerate(axes_flat[:-1]):
#     n_cell = i + 1

#     rmsds = np.load(f"data/rmsds/rmsds_cell{n_cell}.npy")
#     ax.plot(rmsds[-1])

#     # rmsds_ref = calc_rmsds(f"data/trajs/jan/traj_cell{n_cell}.gsd")

#     # ax.plot(rmsds_ref, "C1-")

#     ax.set_xlabel(f"Configurations Cell {n_cell}")
#     ax.set_ylabel("RMSD")
# # y_max = np.max([ax.get_ylim()[1] for ax in axes_flat])

# fig.show()

# # for ax in axes_flat:
# #     ax.set_ylim(0, y_max)

# %% rmsd grouping

# print("grouping\n")

# # n_cell = 5

# # CUTOFF_RMSD = 4


# # rmsds = np.load(f"data/rmsds/rmsds_cell{n_cell}.npy")
# # idxs = set(range(105))

# # groups = []

# # while idxs:
# #     idx = list(idxs)[0]

# #     group = list(idx + (rmsds[idx, idx:] < CUTOFF_RMSD).nonzero()[0])

# #     groups += [group]
# #     idxs -= set(group)
# # for group in groups:
# #     print(group)

# n_cell = 5

# CUTOFF_RMSD = 3

# print(f"grouping for cell {n_cell}")

# print()


# rmsds = np.load(f"data/rmsds/rmsds_cell{n_cell}.npy")
# idxs = set(range(105))

# groups = []

# while idxs:
#     idx = list(idxs)[0]

#     group = []

#     new_group = {idx}

#     while new_group:
#         group += list(new_group)

#         for idx2 in group:
#             new_group |= set((rmsds[idx2, :] < CUTOFF_RMSD).nonzero()[0])
#         new_group -= set(group)
#     groups += [sorted(group)]
#     idxs -= set(group)
# for group in groups:
#     print(group)
# for i in range(len(groups)):
#     for j in range(i):
#         assert not (set(groups[i]) & set(groups[j])), f"Intersection between groups {i} and {j}!"
# %% Plot RMSD grouping

# rmsds = np.load(f"data/rmsds/rmsds_cell{n_cell}.npy")

# # rmsds_last = rmsd.cell_file_rmsds(f"data/trajs/traj_cell5.gsd", ref_idx=-1)

# # rmsds_tenth = rmsd.cell_file_rmsds(f"data/trajs/traj_cell5.gsd", ref_idx=10)

# # rmsds_last = rmsds[-1, :].squeeze()
# # rmsds_tenth = rmsds[10, :].squeeze()

# # fig, (ax0, ax1) = plt.subplots(1, 2, figsize=(14,6))

# # ax0.plot(rmsds_last, "C0-")

# # ax0.set_xlabel(f"Configurations Cell {n_cell}")
# # ax0.set_ylabel("RMSD")


# # ax1.plot(rmsds_tenth, "C0-")

# # ax1.set_xlabel(f"Configurations Cell {n_cell}")
# # ax1.set_ylabel("RMSD")

# fig, ax = new_fig()

# # ax.plot(rmsds_last, "C0-")

# others = []
# del_idxs = []
# for idx, group in enumerate(groups):
#     if len(group) < 5:
#         others += group
#         del_idxs += [idx]
# for del_idx in reversed(del_idxs):
#     del groups[del_idx]
# others_idx = len(groups)

# make_labels = True

# for idx, group in enumerate(groups):
#     # if len(group) < 5:
#     #     continue
#     ref_idx = group[len(group) // 2]  # int(np.median(group))
#     assert ref_idx in group, f"Ref_idx {ref_idx} not in group {idx}!"
#     ax.plot(rmsds[ref_idx, :].squeeze(), f"C{idx}-")  # , label=f"Group {idx} as reference", zorder=-1
#     for idx2, group2 in enumerate(groups):
#         ax.plot(
#             group2,
#             rmsds[ref_idx, :].squeeze()[group2],
#             f"C{idx2}o",
#             markersize=3,
#             label=(f"Group {idx2}" if make_labels else None),
#         )
#     ax.plot(
#         others,
#         rmsds[ref_idx, :].squeeze()[others],
#         f"C{others_idx}o",
#         markersize=3,
#         label=("Others" if make_labels else None),
#     )
#     # Hack to generate labels only once
#     make_labels = False
# # ax.plot(rmsds_last, "C0-", label="Relative to 104")
# # ax.plot(rmsds_tenth, "C1-", label="Relative to 10")

# # ax.plot(groups[0], rmsds_last[groups[0]], "C0o")
# # ax.plot(groups[1], rmsds_last[groups[1]], "C1o")

# # for idx in range(len(groups)):
# #     ax.plot(groups[idx], rmsds_last[groups[idx]], f"C{idx}o", markersize=3, label=f"Group {idx}")

# ax.set_xlabel(f"Configurations Cell {n_cell}")
# ax.set_ylabel("RMSD")

# # ax.legend(loc=(0.6, 0.75))
# ax.legend()

# %% Cell 5 both Configurations

# print("RMSD of both configurations in cell 5")
# print()

# n_cell = 5

# f = gsd.hoomd.open(f"data/cell{n_cell}.gsd")

# pos1 = f[0].particles.position
# pos2 = f[1].particles.position

# r = rmsd.rmsd(pos1, pos2)

# print(r)

# %% Cross-cell RMSDs

print("RMSD of last configuration of each cell")

print()

pos = np.empty((8, NUM_PARTICLES, 3))

for n_cell in range(1, 9):
    traj = gsd.hoomd.open(f"data/trajs/traj_cell{n_cell}.gsd")

    pos[n_cell - 1, :, :] = traj[-1].particles.position
rmsds = np.empty((8, 8))

for i in range(8):
    for j in range(i + 1):
        r = rmsd.rmsd(pos[i, :, :], pos[j, :, :]) if i != j else 0

        rmsds[i, j] = r
        rmsds[j, i] = r
print(rmsds)
