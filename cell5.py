import gsd.hoomd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import tools.rmsd as rmsd
from tools.mg_plot import new_fig, set_styling

# %% Potential Energies

df_log1 = pd.read_csv(f"data/logs/log_cell5.log", sep="\t")
df_log2 = pd.read_csv(f"data/logs/log_cell5_1.log", sep="\t")
df_log3 = pd.read_csv(f"data/logs/log_cell5_2.log", sep="\t")

energies1 = df_log1["potential_energy"]
energies2 = df_log2["potential_energy"]
energies3 = df_log3["potential_energy"]

fig, ax = new_fig()

ax.plot(range(0, 105), energies1, "C0.", label="First simulation")
ax.plot(range(0, 105), energies2, "C1.", label="Second simulation")
ax.plot(range(0, 105), energies3, "C2.", label="Third simulation")

ax.set_xlabel("frame")
ax.set_ylabel("potential energy")

# ax.xaxis.set_major_formatter(FormatStrFormatter('% 1.2f'))
# ax.ticklabel_format(axis="y", style="plain")

ax.legend(loc=(0.6, 0.77))
set_styling(ax)  # , y_loc=(2e5, 1e5)

# %% RMSDs

# rmsd1 = rmsd.cell_file_rmsds("data/trajs/traj_cell5.gsd")
# rmsd2 = rmsd.cell_file_rmsds("data/trajs/traj_cell5_1.gsd")
# rmsd3 = rmsd.cell_file_rmsds("data/trajs/traj_cell5_2.gsd")

# fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, figsize=(18, 6))

# # ax1 = plt.subplot(121)

# # ax2 = plt.subplot(122, sharey=ax1)

# ax1.plot(rmsd1)

# ax2.plot(rmsd2)

# ax3.plot(rmsd3)

# ax1.set_ylabel("RMSD")

# ax1.set_xlabel(r"Configuration first simulation")
# ax2.set_xlabel(r"Configuration second simulation")
# ax3.set_xlabel(r"Configuration third simulation")

# set_styling([ax1, ax2, ax3])


f1 = gsd.hoomd.open("data/trajs/traj_cell5.gsd")
f2 = gsd.hoomd.open("data/trajs/traj_cell5_1.gsd")
f3 = gsd.hoomd.open("data/trajs/traj_cell5_2.gsd")

pos1 = [snap.particles.position for snap in f1]
pos2 = [snap.particles.position for snap in f2]
pos3 = [snap.particles.position for snap in f3]

pos_all = np.stack(pos1 + pos2 + pos3, axis=0)

rmsds = np.empty(pos_all.shape[0])

for i in range(pos_all.shape[0]):
    rmsds[i] = rmsd.rmsd(pos_all[i, :, :], pos_all[104, :, :])

rmsds[104] = 0


fig, ax = new_fig()


ax.plot(range(105), rmsds[:105], "C0-", label="Simulation 1")
ax.plot(range(105, 210), rmsds[105:210], "C1-", label="Simulation 2")
ax.plot(range(210, 315), rmsds[210:315], "C2-", label="Simulation 3")

ax.legend()

ax.set_xlabel("frame")
ax.set_ylabel("RMSD")

set_styling(ax)
