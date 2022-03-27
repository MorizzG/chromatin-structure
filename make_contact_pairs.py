# -*- coding: utf-8 -*-
"""
Created on Sun Jan 30 17:31:30 2022.

@author: mg
"""

# %% Imports

import numpy as np
import pandas as pd

# %% Global Functions


def remove_neighbour_contacts(df_cp):
    """
    Remove neighbouring and self-interacting contacts pairs inplace.

    Parameters
    ----------
    df_contact_pairs : pandas.DataFrame
        data frame containing contact pair information

    """
    df_cp.drop(
        df_cp[
            (df_cp["chr_A"] == df_cp["chr_B"])
            & (
                (df_cp["ind_A"] == df_cp["ind_B"] - 1)
                | (df_cp["ind_A"] == df_cp["ind_B"])
                | (df_cp["ind_A"] == df_cp["ind_B"] + 1)
            )
        ].index,
        inplace=True,
    )


def min_chr_ind(df_cps, n_chr):
    """
    Return smallest ind in a chromosome over all cells.

    Parameters
    ----------
    df : pandas.DataFrame
        Contact data
    ch : int
        number of chromosome

    Returns
    -------
    int
        minimum index

    """
    # return np.min(np.min(df[df["chr_A"] == ch]["ind_A"]), np.min(df[df["chr_B"] == ch]["ind_B"]))

    return np.min(
        [
            (np.min(df_cp[df_cp["chr_A"] == n_chr]["ind_A"]), np.min(df_cp[df_cp["chr_B"] == n_chr]["ind_B"]),)
            for df_cp in df_cps
        ]
    )


def max_chr_ind(df_cps, n_chr):
    """
    Return largest ind in a chromosome over all cells.

    Parameters
    ----------
    df : pandas.DataFrame
        Contact data
    ch : int
        number of chromosome

    Returns
    -------
    int
        maximum index

    """
    # return np.min(np.min(df[df["chr_A"] == ch]["ind_A"]), np.min(df[df["chr_B"] == ch]["ind_B"]))

    return np.max(
        [
            (np.max(df_cp[df_cp["chr_A"] == n_chr]["ind_A"]), np.max(df_cp[df_cp["chr_B"] == n_chr]["ind_B"]),)
            for df_cp in df_cps
        ]
    )


# %% Load CP files and change chromosome identifier from "chrN" to N as number

df_cps = []

for n_cell in range(1, 9):
    df_cp = pd.read_csv(f"data/contact_pairs_raw/Cell{n_cell}_contact_pairs.txt", sep="\t")

    df_cp["chr_A"] = df_cp["chr_A"].apply(lambda s: "chr20" if s == "chrX" else s)
    df_cp["chr_A"] = df_cp["chr_A"].apply(lambda s: int(s[3:]))

    df_cp["chr_B"] = df_cp["chr_B"].apply(lambda s: "chr20" if s == "chrX" else s)
    df_cp["chr_B"] = df_cp["chr_B"].apply(lambda s: int(s[3:]))

    df_cps += [df_cp]
# %% Bin positions with size 100,000 and remove neighbouring contacts

for df_cp in df_cps:
    df_cp["ind_A"] = df_cp["pos_A"] // 100000
    df_cp["ind_B"] = df_cp["pos_B"] // 100000

    remove_neighbour_contacts(df_cp)
# %% Calculate minimum and maximum index and length for each chromosome

chr_mins = pd.Series([min_chr_ind(df_cps, n_chr) for n_chr in range(1, 21)], index=range(1, 21))

chr_maxs = pd.Series([max_chr_ind(df_cps, n_chr) for n_chr in range(1, 21)], index=range(1, 21))

chr_lens = pd.Series([chr_maxs[n_chr] - chr_mins[n_chr] + 1 for n_chr in range(1, 21)], index=range(1, 21),)


chr_lens.to_pickle("data/chromosome_lengths.pkl")

chr_cum_lens = np.cumsum(chr_lens.values)
chr_cum_lens = pd.Series(np.insert(chr_cum_lens, 0, 0)[:-1], index=range(1, 21))


# %% Subtract minimum from ind and add cumsum

# for df_cp in df_cps:
#     df_cp["ind_A"] = df_cp.apply(
#         lambda row: (row["ind_A"] - chr_mins[row["chr_A"]]),
#         axis=1,
#     )
#     df_cp["ind_B"] = df_cp.apply(
#         lambda row: (row["ind_B"] - chr_mins[row["chr_B"]]),
#         axis=1,
#     )

for df_cp in df_cps:
    for n_chr in range(1, 21):
        df_cp.loc[df_cp["chr_A"] == n_chr, "ind_A"] = (
            df_cp[df_cp["chr_A"] == n_chr]["ind_A"] - chr_mins[n_chr] + chr_cum_lens[n_chr]
        )
        df_cp.loc[df_cp["chr_B"] == n_chr, "ind_B"] = (
            df_cp[df_cp["chr_B"] == n_chr]["ind_B"] - chr_mins[n_chr] + chr_cum_lens[n_chr]
        )
    df_cp.sort_values(["chr_A", "chr_B", "ind_A", "ind_B"], inplace=True)
# %% Write contact pairs to pickle

for n_cell in range(1, 9):
    df_cps[n_cell - 1].to_pickle(f"data/contact_pairs/contact_pairs_cell{n_cell}.pkl")
# %% Debug

# df_cp = df_cps[0]

# df_ref = pd.read_pickle("data/contact_pairs_jan/contact_pairs_cell1.pkl")
