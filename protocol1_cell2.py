#!/usr/bin/env python
# coding: utf-8

import datetime

import hoomd.md
import numpy as np
import pandas as pd


class Cube:
    """
    This class is meant to represent a cube, which can be filled with a chosen amount points on random positions
    These particles will later form the initial configuration
    """

    def __init__(self, length):
        self.length = length
        self.random_points = 0

    def generate_points(self):
        self.random_points = np.random.uniform(0, 1, size=(self.length, 3))


def time_now():
    return datetime.datetime.now().strftime("%y-%m-%d-%H-%M-%S")


comment = time_now()  # This will be written as part of the filename and might help do distinguish different results
num_cycles = 105  # The number cycles simulated and therefore the number of conformations obtained

# read a panda dataframe of particles in contact
reading = pd.read_pickle("data/contact_pairs_cell2.pkl")
contact_pairs = reading[["ind_A", "ind_B"]].values
length = pd.read_pickle("data/chromosome_lengths_cell2.pkl")

NContact = contact_pairs.shape[0]  # number of contacts
N_sum = length.sum()  # number of particles
N_sum = int(N_sum)
diff = len(length.keys())  # number of chromosomes, here: 20

# generate the cube for the initial configuration
box = Cube(N_sum)
box.generate_points()

# define the hoomd snapshot here
hoomd.context.initialize("")
s = hoomd.data.make_snapshot(
    N=N_sum,
    box=hoomd.data.boxdim(Lx=50000, Ly=50000, Lz=50000),
    particle_types=["A"],
    bond_types=["backbone", "contact"],
)

s.particles.typeid[:] = 0  # all particles of type A
s.particles.position[:] = box.random_points[:, :]
# Resize the number of HOOMD-Bonds: One from each bead to the next minus 20 since the chromosomes are 20 single chains
# and one for each contact
s.bonds.resize(N_sum - diff + NContact)

# connect all particles of a chromosome to a chain
L = []
n_alreadyadded = 0
for key in sorted(length.keys()):
    a = np.arange(n_alreadyadded, n_alreadyadded + length[key] - 1, dtype="int")
    b = np.arange(n_alreadyadded + 1, n_alreadyadded + length[key], dtype="int")
    n_alreadyadded += length[key]
    temp = np.column_stack((a, b))
    L.append(temp)
L = tuple(L)
final = np.vstack(L)
print(N_sum, " ", final.shape[0], " ", len(length.keys()))
s.bonds.group[: N_sum - diff] = final
s.bonds.group[N_sum - diff :] = contact_pairs  # connect all particles in contact
s.bonds.typeid[: N_sum - diff] = 0  # Set id for chain bonds (bonds)
s.bonds.typeid[N_sum - diff :] = 1  # Set id for contact bonds (contacts)

# setting potentials and their coefficients
system = hoomd.init.read_snapshot(s)  # create a dynamic system from the prepared snapshot
nl = hoomd.md.nlist.tree()  # Set Verlet-list style
# Exclude particles already connected by a bond or contact from the pair interaction defined below
nl.reset_exclusions(exclusions=["bond"])

# define excluded volume interaction
gauss = hoomd.md.pair.gauss(r_cut=0.5, nlist=nl)
gauss.pair_coeff.set("A", "A", epsilon=100.0, sigma=0.5)
gauss.disable(log=True)

# setting the bond potential to harmonic bonds
harmonic = hoomd.md.bond.harmonic()
harmonic.bond_coeff.set("backbone", k=2000.0, r0=1.0)
harmonic.bond_coeff.set("contact", k=2000.0, r0=1.5)

# Chose Langevin thermostat and integrator parameters
all_ = hoomd.group.all()
hoomd.md.integrate.mode_standard(dt=0.001)
seed_langevin = np.random.randint(0, 100000)
integrator = hoomd.md.integrate.langevin(group=all_, kT=1.0, seed=seed_langevin)

# Write out the potential energy in a log and the structures to a gsd after every cycle
per = int(18e4)
hoomd.analyze.log(
    filename="log_all_" + comment + ".log",
    quantities=["potential_energy", "temperature"],
    period=per,
    overwrite=True,
    phase=per - 1,
)
hoomd.dump.gsd(
    "traj_all_" + comment + ".gsd",
    period=per,
    group=all_,
    overwrite=True,
    dynamic=["property", "attribute"],
    phase=per - 1,
)

# start of the simulation
for i in range(1, 1 + num_cycles):
    print(120 * "-")
    print("Step: %i/%i" % (int(i), int(num_cycles)))
    hoomd.run(0.8e5)
    gauss.enable()
    gauss.pair_coeff.set("A", "A", epsilon=100.0, sigma=0.1, r_cut=0.4)
    hoomd.run(0.5e5)
    gauss.pair_coeff.set("A", "A", epsilon=100.0, sigma=1, r_cut=3.5)
    hoomd.run(0.5e5)
    gauss.disable()
