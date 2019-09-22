#!/usr/bin/env python


"""
Script to compare PD190 results obtained in different flowcells
"""

import sys
import os
import string


import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


sns.set(context="paper", style="ticks", palette="colorblind", color_codes=True)
matplotlib.rcParams["svg.fonttype"] = "none"
# Don't use LaTeX for rendering
matplotlib.rcParams["text.usetex"] = False


def main():
    header = "read\tunique_umis\tumi\tgene\thuman\tmouse\ttotal\tmax\tratio\t"
    header += "sp_ratio\tdoublet\tunique_fraction\thuman_norm\ttotal_norm\tmax_norm\t"
    header += "ratio_norm\tsp_ratio_norm\tdoublet_norm\tplate_well"

    exp = {
        "SP": "PD190_humanmouse_sp",
        "S1": "PD190_humanmouse_2"}

    wells = list()
    for _s in string.ascii_uppercase[:8]:
        for _n in [str(x).zfill(2) for x in range(1, 13)]:
            wells.append(_s + _n)

    res = dict()
    for i, (label, _e) in enumerate(exp.items()):
        df = list()
        for well in wells:
            _f = os.path.join(
                "data", _e, f"PD190_humanmouse_{well}",
                f"PD190_humanmouse_{well}.metrics.csv.gz")
            df.append(pd.read_csv(_f, header=None, names=header.split("\t")))
        res[label] = pd.concat(df).sort_values('read')

    fig, axis = plt.subplots(1, 2, figsize=(3 * 2, 3), sharex=True, sharey=False)
    for i, (label, _) in enumerate(exp.items()):
        df = res[label].sort_values('read').tail(1000000).sample(frac=0.1)
        axis[i].scatter(
            df['read'], df['unique_fraction'],
            alpha=0.05, s=2,
            c=df['sp_ratio_norm'], cmap="RdBu_r", vmin=0, vmax=1,
            rasterized=True)
        axis[i].set_xscale("log")
        axis[i].set_ylim((-0.1, 1.1))
        axis[i].set_title(label)
        axis[i].set_xlabel("Reads per cell")
        if i == 0:
            axis[i].set_ylabel("Unique fraction")
        for _x in range(1, 5):
            axis[i].axvline(10 ** _x, linestyle="--", color="grey", linewidth=0.25)
        for _x in range(1, 11):
            axis[i].axhline(0.1 * _x, linestyle="--", color="grey", linewidth=0.25)
    kwargs = {"dpi": 300, "bbox_inches": "tight"}
    fig.savefig(os.path.join("results", "PD190_humanmouse.SP_vs_S1.svg"), **kwargs)
    fig.savefig(os.path.join("results", "PD190_humanmouse.SP_vs_S1.pdf"), **kwargs)


if __name__ == "__main__":
    sys.exit(main())
