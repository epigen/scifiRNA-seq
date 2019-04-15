import sys
import os
import time
import pickle

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import scipy


def main():
    global output_prefix
    output_prefix = os.path.join("results", "joint_analysis.")

    # sample annotations
    sample1 = {
        "sample_name": "sci-RNA-seq_SCI020_3uL_reseq_4K",
        "n_parts": 4,
        "n_round1_barcodes": 96,
        "round2_barcodes": "reverse_complement",
        "expected_cell_number": 4000,
        "r1_barcode_annotation": os.path.join("metadata", "sciRNA-seq.SCI017.oligos_2019-02-11.csv"),
        "r1_barcode_orientation": "barcode_sequence",
        "r2_barcode_annotation": os.path.join("metadata", "737K-cratac-v1.reverse_complement.csv"),
        "r2_barcode_orientation": "original"}

    sample2 = {
        "sample_name": "sci-RNA-seq_SCI021_125K",
        "n_parts": 2,
        "n_round1_barcodes": 96,
        "round2_barcodes": "original",
        "expected_cell_number": 125000,
        "r1_barcode_annotation": os.path.join("metadata", "sciRNA-seq.SCI017.oligos_2019-02-11.csv"),
        "r1_barcode_orientation": "barcode_sequence",
        "r2_barcode_annotation": os.path.join("metadata", "737K-cratac-v1.reverse_complement.csv"),
        "r2_barcode_orientation": "original"}
    samples = [sample1, sample2]

    # convenience
    for sample in samples:
        sample["sample_dir"] = os.path.join("data", sample['sample_name'])
        sample["results_dir"] = "results"
        sample["input_prefix"] = os.path.join(sample["sample_dir"], f"{sample['sample_name']}.")
        sample["output_prefix"] = os.path.join(sample["results_dir"], f"{sample['sample_name']}.")
        sample['metrics'] = pickle.load(open(sample['input_prefix'] + f"metrics.pickle", 'rb'))
        for bc in ["r1", "r2"]:
            sample[bc + "_barcode_sequence"] = pd.read_csv(sample[bc + "_barcode_annotation"])[bc + "_barcode_orientation"]

    # # Plot
    plot_metrics_lineplot(samples, tail=int(1e6), suffix="ar_2.0", aspect_ratio=2)
    plot_metrics_lineplot(samples, tail=int(1e6), suffix="ar_1.1", aspect_ratio=1.1)


def plot_metrics_lineplot(samples, keys=['read', 'umi', 'gene'], tail=None, suffix="", aspect_ratio=2):
    def min_max(x):
        return (x - x.min()) / (x.max() - x.min())

    print(f"# {time.asctime()} - Plotting metrics per cell.")
    n = len(keys)

    colors = sns.color_palette("colorblind", len(samples))

    fig, axis = plt.subplots(n, 2, figsize=(2 * 3 * aspect_ratio, n * 3), tight_layout=True)
    for i, metric in enumerate(keys):
        for j, sample in enumerate(samples):
            m = sample["metrics"]
            if tail:
                t = min(tail, m.shape[0])
                m = m.tail(t)
            d = m.loc[:, metric].sort_values()
            # d_norm = d / d.sum()
            d_norm = min_max(d)
            rank = d.rank(ascending=False, method="average", pct=False)
            # rank_norm = rank / rank.shape[0]
            rank_norm = min_max(rank)
            e = sample['expected_cell_number']
            kwargs = {"rasterized": True, "label": sample['sample_name']}

            # first plot the "real" cells
            axis[i, 0].plot(rank[-e:], d[-e:], color=colors[j], **kwargs)
            axis[i, 1].plot(rank_norm[-e:], d_norm[-e:], color=colors[j], **kwargs)

            # plot lines with expected cell number
            axis[i, 0].axvline(e, linestyle="--", color=colors[j])
            axis[i, 1].axvline((e - rank.min()) / (rank.max() - rank.min()), linestyle="--", color=colors[j])

            # then, the remaining in grey
            axis[i, 0].plot(rank[:-e], d[:-e], color="grey", **kwargs)
            axis[i, 1].plot(rank_norm[:-e], d_norm[:-e], color="grey", **kwargs)
        for ax in axis[i, :]:
            ax.loglog()
            ax.set_xlabel("Barcodes")
            ax.set_ylabel(metric.capitalize() + "s")
            ax.legend()
    fig.savefig(
        output_prefix + f"metrics_per_cell.lineplot.{suffix}.svg"
        .replace("..", "."),
        dpi=300, bbox_inches="tight")


def plot_metrics_lineplot(samples, keys=['read', 'umi', 'gene'], tail=None, suffix="", aspect_ratio=2):
    def min_max(x):
        return (x - x.min()) / (x.max() - x.min())

    print(f"# {time.asctime()} - Plotting metrics per cell.")
    n = len(keys)

    colors = sns.color_palette("colorblind", len(samples))

    fig, axis = plt.subplots(n, 1, figsize=(1 * 3 * aspect_ratio, n * 3), tight_layout=True, squeeze=False)
    for i, metric in enumerate(keys):
        for j, sample in enumerate(samples):
            m = sample["metrics"]
            if tail:
                t = min(tail, m.shape[0])
                m = m.tail(t)
            d = m.loc[:, metric].sort_values()
            # d_norm = d / d.sum()
            d_norm = min_max(d)
            rank = d.rank(ascending=False, method="average", pct=False)
            # rank_norm = rank / rank.shape[0]
            rank_norm = min_max(rank)
            e = sample['expected_cell_number']
            kwargs = {"rasterized": True, "label": sample['sample_name']}

            # first plot the "real" cells
            axis[i, 0].plot(rank[-e:], d[-e:], color=colors[j], **kwargs)  # real cells
            axis[i, 0].plot(rank[:-e], d[:-e], color="grey", **kwargs)  # empty barcodes
            ax2 = axis[i, 0].twinx()
            ax2.plot(rank_norm[-e:], d_norm[-e:], color=colors[j], visible=False, **kwargs)
            ax2.plot(rank_norm[:-e], d_norm[:-e], color="grey", visible=False, **kwargs)

            # plot lines with expected cell number
            # axis[i, 0].axvline(e, linestyle="--", color=colors[j])
            # ax2.axvline((e - rank.min()) / (rank.max() - rank.min()), linestyle="--", color=colors[j])

        for ax in axis[i, :]:
            ax.loglog()
            ax.set_xlabel("Barcodes")
            ax.set_ylabel(metric.capitalize() + "s")
            ax.legend()
    fig.savefig(
        output_prefix + f"metrics_per_cell.lineplot.{suffix}.svg"
        .replace("..", "."),
        dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    sns.set(context="paper", style="ticks", palette="colorblind", color_codes=True)
    matplotlib.rcParams["svg.fonttype"] = "none"
    # Don't use LaTeX for rendering
    matplotlib.rcParams["text.usetex"] = False
    sys.exit(main())
