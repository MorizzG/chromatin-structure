import numpy as np
import pandas as pd


def calc_chrom_lens(df_contact_pairs):

    chr_min_pos = np.empty(20)
    chr_max_pos = np.empty(20)

    for n in range(1, 21):
        chrn_minA = min(df_contact_pairs.loc[df_contact_pairs["chr_A"] == n]["pos_A"])
        chrn_minB = min(df_contact_pairs.loc[df_contact_pairs["chr_B"] == n]["pos_B"])

        chr_min_pos[n - 1] = min(chrn_minA, chrn_minB)

        chrn_maxA = max(df_contact_pairs.loc[df_contact_pairs["chr_A"] == n]["pos_A"])
        chrn_maxB = max(df_contact_pairs.loc[df_contact_pairs["chr_B"] == n]["pos_B"])

        chr_max_pos[n - 1] = max(chrn_maxA, chrn_maxB)

        return chr_max_pos - chr_min_pos


def make_contact_pairs(df_contact_pairs):
    # chromosomes are initially labelled as "chr17" or "chrX" for ch 20
    # we want them labelled by numbers instead: chr17 -> 17, chrX -> 20
    df_contact_pairs["chr_A"] = df_contact_pairs["chr_A"].apply(lambda s: "chr20" if s == "chrX" else s)
    df_contact_pairs["chr_A"] = df_contact_pairs["chr_A"].apply(lambda s: int(s[3:]))

    df_contact_pairs["chr_B"] = df_contact_pairs["chr_B"].apply(lambda s: "chr20" if s == "chrX" else s)
    df_contact_pairs["chr_B"] = df_contact_pairs["chr_B"].apply(lambda s: int(s[3:]))

    chr_min_pos = {}

    # calculate min pos for each chromosome
    for n in range(1, 21):
        chrn_minA = min(df_contact_pairs.loc[df_contact_pairs["chr_A"] == n]["pos_A"])
        chrn_minB = min(df_contact_pairs.loc[df_contact_pairs["chr_B"] == n]["pos_B"])

        chr_min_pos[n] = min(chrn_minA, chrn_minB)

    # lengths = pd.read_pickle("data/chromosome_lengths.pkl")

    lengths = calc_chrom_lens(df_contact_pairs)

    print(lengths)
    return

    lengths_cumsum = np.array(np.cumsum(lengths))
    lengths_cumsum = np.insert(lengths_cumsum, 0, 0)

    # convert position on chromosomes to bead number
    # each bead is a bin of 100,000 chromosome positions
    # subtract the minimin position so our beads start at 0
    # which chromosome a bead belongs to will be later taken care of
    df_contact_pairs["ind_A"] = df_contact_pairs.apply(
        lambda row: (row["pos_A"] - chr_min_pos[row["chr_A"]]) // 100000 + lengths_cumsum[row["chr_A"] - 1], axis=1,
    )
    df_contact_pairs["ind_B"] = df_contact_pairs.apply(
        lambda row: (row["pos_B"] - chr_min_pos[row["chr_B"]]) // 100000 + lengths_cumsum[row["chr_B"] - 1], axis=1,
    )

    # df_contact_pairs["ind_A"] = df_contact_pairs.apply(
    #     lambda row: (row["pos_A"] - 3000000) // 100000 + lengths_cumsum[row["chr_A"] - 1],
    #     axis = 1,
    # )
    # df_contact_pairs["ind_B"] = df_contact_pairs.apply(
    #     lambda row: (row["pos_B"] - 3000000) // 100000 + lengths_cumsum[row["chr_B"] - 1],
    #     axis = 1,
    # )

    # drop self-interacting and neighbouring terms
    df_contact_pairs = df_contact_pairs.loc[
        (df_contact_pairs["ind_A"] != df_contact_pairs["ind_B"])
        & (df_contact_pairs["ind_A"] != df_contact_pairs["ind_B"] + 1)
        & (df_contact_pairs["ind_A"] != df_contact_pairs["ind_B"] - 1)
    ]

    return df_contact_pairs[["ind_A", "ind_B"]].values
