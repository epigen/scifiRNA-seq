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
    global sample_name
    global sample_dir
    global results_dir
    global input_prefix
    global output_prefix
    global expected_cell_number
    global n_round1_barcodes

    # barcode annotations
    annotation_file = os.path.join("metadata", "sciRNA-seq.SCI022.oligos_2019-04-14.csv")
    annotation = pd.read_csv(annotation_file)
    annotation.loc[:, "sample_type"] = pd.np.nan
    annotation.loc[annotation['sample'].str.startswith("PBMC"), "sample_type"] = "PBMC"
    annotation.loc[annotation['sample'].str.startswith("Tcell"), "sample_type"] = "TCell"
    annotation.loc[annotation['sample'].str.startswith("3T3"), "sample_type"] = "Mixture"

    droplet_barcode_file = os.path.join("metadata", "737K-cratac-v1.reverse_complement.csv")
    droplets = pd.read_csv(droplet_barcode_file)

    # sample annotations

    sample_name = "SCI022_TCell"
    n_lanes = 2
    n_parts = 4
    n_round1_barcodes = 384
    round2_barcodes = "original"
    expected_cell_number = 125000

    sample_name = "SCI022_PBMC"
    n_lanes = 2
    n_parts = 4
    n_round1_barcodes = 384
    round2_barcodes = "original"
    expected_cell_number = 125000

    # convenience
    sample_dir = os.path.join("data", sample_name)
    results_dir = "results"
    input_prefix = os.path.join(sample_dir, f"{sample_name}.")
    output_prefix = os.path.join(results_dir, f"{sample_name}.")

    # read text files
    df = load_data(sample_name, n_lanes=n_lanes, n_parts=n_parts)

    print(f"# {time.asctime()} - Writing data to pickle.")
    # pickle.dump(df, open(input_prefix + "all.pickle", 'wb'), protocol=-1)
    # df = pickle.load(open(input_prefix + "all.pickle", 'rb'))

    # Gather metrics per cell
    r1_annotation = annotation.set_index("barcode_sequence")['sample_type']
    r1_annotation.index.name = "r1"
    metrics = gather_stats_per_cell(df, r1_annotation=r1_annotation)

    # # Plot
    # # # Loglog line plot of rank vs abundance
    t = int(expected_cell_number * 10)
    plot_metrics_lineplot(metrics, tail=t)
    plot_metrics_lineplot(metrics, tail=None, by_group="sample_type")
    # # # Histogram
    t = int(expected_cell_number * 10)
    plot_metrics_distplot(metrics, tail=t)

    # # # Efficiency plot (UMI and gene yield in relation to reads per cell)
    t = int(expected_cell_number * 5)
    plot_efficiency(metrics, tail=t)
    plot_efficiency(metrics, tail=t, by_group="sample_type")
    plot_efficiency(metrics, tail=t, colour_by="unique_fraction")

    # # # Species mixing - only for twice the number of expected cells
    t = int(expected_cell_number * 1.5)
    plot_species_mixing(metrics, tail=t, norm=False)
    plot_species_mixing(metrics, tail=t, norm=False, by_group="sample_type")
    # plot_species_mixing(metrics, tail=t, norm=True)
    # plot_species_mixing(metrics, tail=t, norm=True, by_group="sample_type")

    # Well and droplet inspection
    well_metrics = gather_stats_per_well(metrics, seq_content=True)
    well_metrics.to_csv(output_prefix + "well_metrics.csv.gz")
    plot_well_stats(well_metrics, tail=None, suffix="")

    #

    # Now do the same only for exactly matching barcodes
    metrics_filtered = get_exact_matches(
        metrics,
        r1_whitelist=annotation['barcode_sequence'],
        r2_whitelist=droplets[round2_barcodes],
        save_intermediate=True, plot=True)

    t = expected_cell_number * 10
    plot_metrics_lineplot(metrics_filtered, tail=t, suffix="exact_match")
    plot_metrics_lineplot(metrics_filtered, tail=None, suffix="exact_match", by_group="sample_type")

    t = expected_cell_number * 10
    plot_metrics_distplot(metrics_filtered, tail=t, suffix="exact_match")
    plot_metrics_distplot(metrics_filtered, tail=t, suffix="exact_match", by_group="sample_type")

    # # # Efficiency plot (UMI and gene yield in relation to reads per cell)
    t = expected_cell_number * 5
    plot_efficiency(metrics_filtered, tail=t, suffix="exact_match")
    plot_efficiency(metrics_filtered, tail=t, suffix="exact_match", colour_by="unique_fraction")
    plot_efficiency(metrics_filtered, tail=t, suffix="exact_match", by_group="sample_type")


    # # # Species mixing - only for twice the number of expected cells
    t = int(expected_cell_number * 1.5)
    plot_species_mixing(metrics_filtered, tail=t, norm=False, suffix="exact_match")
    # plot_species_mixing(metrics_filtered, tail=t, norm=False, suffix="exact_match", by_group="sample_type")
    # plot_species_mixing(metrics_filtered, tail=t, norm=True, suffix="exact_match")
    # plot_species_mixing(metrics_filtered, tail=t, norm=True, suffix="exact_match", by_group="sample_type")

    # Well and droplet inspection
    well_metrics_filtered = gather_stats_per_well(metrics_filtered, seq_content=True, save_intermediate=True)
    well_metrics_filtered.to_csv(output_prefix + "well_metrics_filtered.csv.gz")
    plot_well_stats(well_metrics_filtered, tail=None, suffix="exact_match")

    #
    metrics_r2 = gather_stats_per_cell_as_droplet(
        df,
        cell_barcodes=['r2'],
        doublet_threshold=0.85,
        save_intermediate=True,
        norm_species=True,
        r1_annotation=None,
        suffix="only_r2")
    metrics_r2_filtered = get_exact_matches_droplet(metrics_r2, droplets[round2_barcodes])

    # # assign singlet or doublet across range of required purity
    for (o, d), suffix in [
            ((metrics, metrics_r2), ""),
            ((metrics_filtered, metrics_r2_filtered), "exact_match")
    ]:
        c = o.query("total > 100")
        d = d.query("total > 100")
        purity = dict()
        droplet_purity = dict()
        for t in np.arange(0.4, 1.00, 0.001):
            purity[t] = (c.loc[:, 'ratio'] < t).sum() / c.shape[0]
            droplet_purity[t] = (d.loc[:, 'ratio'] < t).sum() / d.shape[0]
            # print(f"At 't': {t:.2f}, doublet rate: {purity[t]:.3f}")

        purity = pd.Series(purity, name="doublet_rate")
        purity.index.name = "purity"
        purity[1.0] = 1.0
        droplet_purity = pd.Series(droplet_purity, name="doublet_rate")
        droplet_purity.index.name = "purity"
        droplet_purity[1.0] = 1.0

        fig, axis = plt.subplots(1, figsize=(3, 3), tight_layout=True)
        for x in np.linspace(0.5, 1.0, 6):
            axis.axvline(x, linestyle="--", color="lightgrey")
        for y in np.linspace(0.0, 1.0, 11):
            axis.axhline(y, linestyle="--", color="lightgrey")
        axis.set_xlim((0.4, 1.1))
        axis.plot(purity.index, purity, label="scifi-RNA-seq")
        axis.plot(droplet_purity.index, droplet_purity, label="10X only")
        axis.legend()
        axis.set_xlabel("Purity required\n(fraction of transcriptome in genome)")
        axis.set_ylabel("Doublet rate")
        fig.savefig(
            output_prefix + f".umis_per_cell.purity_range.comparison_to_droplet.{suffix}._r2_only.lineplot.svg",
            dpi=300, bbox_inches="tight")

    # Plot species mixing compared to 10X only barcodes
    c = metrics_filtered.query("total > 100").tail(2 * expected_cell_number)
    d = metrics_r2_filtered.query("total > 100").tail(2 * expected_cell_number)
    plot_comparison_to_10x(c, d, suffix="_r2_only")

    #

    # Investigate cells per droplet
    # # # Count
    cells_per_droplet = metrics_filtered.query("umi > 250").reset_index().groupby(['r2'])['r1'].nunique()
    cells_per_droplet.name = 'cells_per_droplet'

    cells_per_droplet_counts = cells_per_droplet.value_counts().sort_index()

    # # # Try to estimate zero
    zero = metrics_filtered.loc[(metrics_filtered['umi'] < 100) & (metrics_filtered['umi'] >= 10)].index.get_level_values("r2").unique().shape[0]
    cells_per_droplet_counts[0] = zero
    cells_per_droplet_counts = cells_per_droplet_counts.sort_index()
    # add pseudocounts to observed
    cells_per_droplet = cells_per_droplet.append(pd.Series(np.zeros(zero, dtype=int)))

    to_pickle(cells_per_droplet_counts, "cells_per_droplet_counts")
    cells_per_droplet_counts.to_csv(output_prefix + "cells_per_droplet_counts.csv", header=True, index=True)

    # # Investigate statistical properties
    cells_per_droplet_stats(cells_per_droplet)

    #
    # # Investigate dependency of cells per droplet on transcriptome size
    metrics_tmp = metrics_filtered.join(cells_per_droplet)

    metrics_tmp2 = metrics_tmp.tail(10000)

    fig, axis = plt.subplots(2, 3, figsize=(3 * 3, 2 * 3), tight_layout=True)
    sns.distplot(metrics_tmp2['cells_per_droplet'], bins=40, kde=False, ax=axis[0, 0])
    axis[0, 1].scatter(
        metrics_tmp2['cells_per_droplet'], metrics_tmp2['umi'],
        rasterized=True, alpha=0.1, s=2)
    sns.boxplot(data=metrics_tmp2, x="cells_per_droplet", y='umi', ax=axis[0, 2], whis=10000)

    stats2 = metrics_tmp.query("umi > 500")

    sns.distplot(stats2['cells_per_droplet'], bins=40, kde=False, ax=axis[1, 0])
    axis[1, 1].scatter(
        stats2['cells_per_droplet'], stats2['umi'],
        rasterized=True, alpha=0.1, s=2)
    sns.boxplot(data=stats2, x="cells_per_droplet", y='umi', ax=axis[1, 2], whis=10000)
    for ax in axis.flatten():
        ax.set_yscale("log")
    axis[0, 1].set_xscale("log")
    axis[1, 1].set_xscale("log")
    fig.savefig(
        output_prefix + f"cells_per_droplet_vs_transcriptome_size.svg",
        dpi=300, bbox_inches="tight")


def load_data(sample_name, n_lanes=2, n_parts=4, **kwargs):
    print(f"# {time.asctime()} - Reading files.")
    pieces = list()
    for lane in range(1, n_lanes + 1):
        for part in range(1, n_parts + 1):
            print(f"# {time.asctime()} - Reading lane {lane}, part {part}.")
            input_file = os.path.join(f"data/{sample_name}/{sample_name}.{lane}.{part}.STAR.filtered2.csv.gz")

            d = pd.read_csv(
                input_file,
                header=None, sep=',', error_bad_lines=False, engine='c',
                names=['read', 'r1', 'r2', 'umi', 'gene', 'pos'], compression='gzip', **kwargs)
            pieces.append(d)
            print(f"# {time.asctime()} - Done with lane {lane}, part {part}. {d.shape[0]} lines.")

    print(f"# {time.asctime()} - Concatenating parts.")
    return pd.concat(pieces)


def gather_stats_per_cell(
        df,
        cell_barcodes=['r1', 'r2'],
        doublet_threshold=0.85,
        save_intermediate=True,
        norm_species=True,
        r1_annotation=None,
        suffix=""):
    print(f"# {time.asctime()} - Gathering metrics per cell.")

    # performance metrics
    # # number of unique reads per cell
    reads_per_cell = df.groupby(['r1', 'r2'], sort=False)['read'].nunique().sort_values()
    if save_intermediate:
        to_pickle(reads_per_cell, 'reads_per_cell' + suffix)

    # # mapping rate per cell
    # TODO: add mapping rate per cell (needs additional file)

    # # duplication per cell
    reads_per_umi = df.groupby(cell_barcodes + ['gene', 'pos'], sort=False)['umi'].size()
    reads_per_umi = reads_per_umi.reset_index(level=['gene', 'pos'], drop=True)
    if save_intermediate:
        to_pickle(reads_per_umi, 'reads_per_umi' + suffix)

    unique = (reads_per_umi == 1)
    unique_per_cell = unique.groupby(level=cell_barcodes).sum()
    unique_per_cell.name = "unique_umis"
    if save_intermediate:
        to_pickle(unique_per_cell, 'unique_per_cell' + suffix)

    # # UMI count per cell
    umi_counts = df.groupby(cell_barcodes + ['gene', 'pos'], sort=False)['umi'].nunique()
    if save_intermediate:
        to_pickle(umi_counts, 'umi_counts' + suffix)
    umis_per_cell = umi_counts.groupby(cell_barcodes, sort=False).sum().sort_values()
    if save_intermediate:
        to_pickle(umis_per_cell, 'umis_per_cell' + suffix)

    # # Genes per cell
    genes_per_cell = umi_counts.reset_index(level='gene').groupby(cell_barcodes, sort=False)['gene'].nunique().sort_values()
    if save_intermediate:
        to_pickle(genes_per_cell, 'genes_per_cell' + suffix)

    # species mixing
    print(f"# {time.asctime()} - Gathering species-specific metrics per cell.")
    umi_counts2 = umi_counts.reset_index()
    umi_counts2 = (
        umi_counts2
        .assign(
            species=(
                umi_counts2['gene'].str.startswith("ENSG")
                .replace(True, "human").replace(False, "mouse"))))
    species_counts = umi_counts2.groupby(cell_barcodes + ['species'])['umi'].sum()

    spc = species_counts.reset_index().pivot_table(
        index=cell_barcodes, columns='species', values='umi', fill_value=0)
    spc += 1
    spc = spc.assign(total=spc.sum(1), max=spc.max(1))
    spc = (
        spc.assign(
            ratio=spc["max"] / spc['total'],
            sp_ratio=spc["human"] / spc['total'])
        .sort_values("total"))
    spc = spc.assign(doublet=(spc['ratio'] < doublet_threshold).astype(int).replace(0, -1))
    if save_intermediate:
        to_pickle(spc, "spc" + suffix)

    print(f"# {time.asctime()} - Joining all metrics.")
    # TODO: speed up by assigning to column using equallly sorted indexes
    metrics = (
        reads_per_cell.to_frame()
        .join(unique_per_cell)
        .join(umis_per_cell)
        .join(genes_per_cell)
        .join(spc).sort_values("umi"))

    # Calculate unique fraction
    metrics.loc[:, "unique_fraction"] = metrics['unique_umis'] / metrics['umi']

    if r1_annotation is not None:
        print(f"# {time.asctime()} - Adding well annotation to metrics.")
        r1_annotation.index.name = "r1"
        metrics = metrics.join(r1_annotation)
    if save_intermediate:
        to_pickle(metrics, "metrics" + suffix)

    if not norm_species:
        return metrics

    # Add normalized metrics to stats
    metrics.loc[:, 'human_norm'] = metrics['human'] * r[int(expected_cell_number)]
    metrics.loc[:, 'total_norm'] = metrics[['mouse', 'human_norm']].sum(1)
    metrics.loc[:, 'max_norm'] = metrics[['mouse', 'human_norm']].max(1)
    metrics.loc[:, 'ratio_norm'] = metrics['max_norm'] / metrics['total_norm']
    metrics.loc[:, 'sp_ratio_norm'] = metrics['human_norm'] / metrics['total_norm']
    metrics.loc[:, 'doublet_norm'] = (metrics.loc[:, 'ratio_norm'] < doublet_threshold).astype(int).replace(0, -1)

    if save_intermediate:
        to_pickle(metrics, "metrics" + suffix)

    # Assess species bias
    r = dict()
    for f in [0.1, 0.2, 0.25, 0.5, 0.75] + list(range(1, 1000, 10)) + [1.5]:
        t = metrics.tail(int(expected_cell_number * f))
        r[expected_cell_number * f] = t['mouse'].mean() / t['human'].mean()
    r = pd.Series(r).sort_index()

    # # plot ammount of species-specific bias
    n = 1
    fig, axis = plt.subplots(n, 1, figsize=(3, n * 3), tight_layout=True)
    axis.plot(r.index, r, "o-")
    axis.axvline(expected_cell_number, linestyle="--", color="grey")
    axis.set_xlabel("Top N barcodes taken into account")
    axis.set_ylabel("Ratio of mean UMIs per cell\nbetween mouse and human")
    axis.set_xscale("log")
    fig.savefig(
        output_prefix + f"species_bias.lineplot.svg"
        .replace("..", "."),
        dpi=300, bbox_inches="tight")

    # Plot illustration of normalization procedure
    metrics2 = metrics.tail(int(2 * expected_cell_number))

    for label in ["", ".log"]:
        fig, axis = plt.subplots(2, 2, figsize=(2 * 3, 2 * 3), tight_layout=True)
        kwargs = {"s": 0.1, "alpha": 0.2, "rasterized": True, "cmap": get_custom_cmap()}
        axis[0, 0].scatter(
            metrics2['mouse'],
            metrics2['human'],
            c=metrics2['sp_ratio'], **kwargs)
        axis[0, 1].scatter(
            metrics2['mouse'],
            metrics2['human'],
            c=metrics2['sp_ratio_norm'], **kwargs)
        v = metrics2[['mouse', 'human']].quantile(.999).max()
        v += v * 0.1
        for ax in axis[0, :]:
            if label == "":
                ax.set_xlim((-(v / 30.), v))
                ax.set_ylim((-(v / 30.), v))
            else:
                ax.loglog()
            ax.plot((0, v), (0, v), linestyle="--", color="grey")
            ax.plot((0, v), (0, v * (1 / r[int(expected_cell_number)])), linestyle="--", color="orange")
            ax.set_ylabel("Human (UMIs)")

        axis[1, 0].scatter(
            metrics2['mouse'],
            metrics2['human_norm'],
            c=metrics2['sp_ratio'], **kwargs)
        axis[1, 1].scatter(
            metrics2['mouse'],
            metrics2['human_norm'],
            c=metrics2['sp_ratio_norm'], **kwargs)
        v = metrics2[['mouse', 'human_norm']].quantile(.999).max()
        v += v * 0.1
        for ax in axis[1, :]:
            if label == "":
                ax.set_xlim((-(v / 30.), v))
                ax.set_ylim((-(v / 30.), v))
            else:
                ax.loglog()
            ax.plot((0, v), (0, v), linestyle="--", color="grey")
            ax.plot((0, v), (0, v), linestyle="--", color="orange")
            ax.set_ylabel("Human (norm UMIs)")

        for ax in axis.flatten():
            ax.set_xlabel("Mouse (UMIs)")
        axis[0, 0].set_title("Coloured by ratio")
        axis[0, 1].set_title("Coloured by norm ratio")
        fig.savefig(
            output_prefix + f"species_bias.original_vs_corrected.scatterplot{label}.svg"
            .replace("..", "."),
            dpi=300, bbox_inches="tight")

    return metrics


def gather_stats_per_cell_as_droplet(
        df,
        cell_barcodes=['r2'],
        doublet_threshold=0.85,
        save_intermediate=True,
        norm_species=True,
        r1_annotation=None,
        suffix="_r2_only"):
    print(f"# {time.asctime()} - Gathering metrics per cell.")

    # performance metrics
    # # UMI count per cell
    umi_counts = df.groupby(cell_barcodes + ['gene', 'pos'], sort=False)['umi'].nunique()
    if save_intermediate:
        to_pickle(umi_counts, 'umi_counts' + suffix)
    umis_per_cell = umi_counts.groupby(cell_barcodes, sort=False).sum().sort_values()
    if save_intermediate:
        to_pickle(umis_per_cell, 'umis_per_cell' + suffix)

    # species mixing
    print(f"# {time.asctime()} - Gathering species-specific metrics per cell.")
    umi_counts2 = umi_counts.reset_index()
    umi_counts2 = (
        umi_counts2
        .assign(
            species=(
                umi_counts2['gene'].str.startswith("ENSG")
                .replace(True, "human").replace(False, "mouse"))))
    species_counts = umi_counts2.groupby(cell_barcodes + ['species'])['umi'].sum()

    spc = species_counts.reset_index().pivot_table(
        index=cell_barcodes, columns='species', values='umi', fill_value=0)
    spc += 1
    spc = spc.assign(total=spc.sum(1), max=spc.max(1))
    spc = (
        spc.assign(
            ratio=spc["max"] / spc['total'],
            sp_ratio=spc["human"] / spc['total'])
        .sort_values("total"))
    spc = spc.assign(doublet=(spc['ratio'] < doublet_threshold).astype(int).replace(0, -1))
    if save_intermediate:
        to_pickle(spc, "spc" + suffix)

    print(f"# {time.asctime()} - Joining all metrics.")
    # TODO: speed up by assigning to column using equallly sorted indexes
    metrics_r2 = (
        umis_per_cell.to_frame()
        .join(spc).sort_values("umi"))

    if save_intermediate:
        to_pickle(metrics_r2, "metrics_r2" + suffix)

    if not norm_species:
        return metrics_r2

    # Assess species bias
    r = dict()
    for f in [0.1, 0.2, 0.25, 0.5, 0.75] + list(range(1, 1000, 10)) + [1.5]:
        t = metrics_r2.tail(int(expected_cell_number * f))
        r[expected_cell_number * f] = t['mouse'].mean() / t['human'].mean()
    r = pd.Series(r).sort_index()

    # Add normalized metrics to stats
    metrics_r2.loc[:, 'human_norm'] = metrics_r2['human'] * r[int(expected_cell_number)]
    metrics_r2.loc[:, 'total_norm'] = metrics_r2[['mouse', 'human_norm']].sum(1)
    metrics_r2.loc[:, 'max_norm'] = metrics_r2[['mouse', 'human_norm']].max(1)
    metrics_r2.loc[:, 'ratio_norm'] = metrics_r2['max_norm'] / metrics_r2['total_norm']
    metrics_r2.loc[:, 'sp_ratio_norm'] = metrics_r2['human_norm'] / metrics_r2['total_norm']
    metrics_r2.loc[:, 'doublet_norm'] = (metrics_r2.loc[:, 'ratio_norm'] < doublet_threshold).astype(int).replace(0, -1)

    if save_intermediate:
        to_pickle(metrics_r2, "metrics_r2" + suffix)

    # Plot ammount of species-specific bias
    n = 1
    fig, axis = plt.subplots(n, 1, figsize=(3, n * 3), tight_layout=True)
    axis.plot(r.index, r, "o-")
    axis.axvline(expected_cell_number, linestyle="--", color="grey")
    axis.set_xlabel("Top N barcodes taken into account")
    axis.set_ylabel("Ratio of mean UMIs per cell\nbetween mouse and human")
    axis.set_xscale("log")
    fig.savefig(
        output_prefix + f"species_bias.lineplot.{suffix}.svg"
        .replace("..", "."),
        dpi=300, bbox_inches="tight")

    # Plot illustration of normalization procedure
    metrics2 = metrics_r2.tail(int(10 * expected_cell_number))

    for label in ["", ".log"]:
        fig, axis = plt.subplots(1, 2, figsize=(2 * 3, 1 * 3), tight_layout=True, squeeze=False)
        kwargs = {"s": 0.1, "alpha": 0.2, "rasterized": True, "cmap": get_custom_cmap()}
        axis[0, 0].scatter(metrics2['mouse'], metrics2['human'], c=metrics2['sp_ratio'], **kwargs)
        # axis[0, 1].scatter(metrics2['mouse'], metrics2['human'], c=metrics2['sp_ratio_norm'], **kwargs)
        v = metrics2[['mouse', 'human']].quantile(.999).max()
        v += v * 0.1
        for ax in axis[0, :]:
            if label == "":
                ax.set_xlim((-(v / 30.), v))
                ax.set_ylim((-(v / 30.), v))
            else:
                ax.loglog()
            ax.plot((0, v), (0, v), linestyle="--", color="grey")
            ax.plot((0, v), (0, v * (1 / r[int(expected_cell_number)])), linestyle="--", color="orange")
            ax.set_ylabel("Human (UMIs)")

        for ax in axis.flatten():
            ax.set_xlabel("Mouse (UMIs)")
        axis[0, 0].set_title("Coloured by ratio")
        axis[0, 1].set_title("Coloured by norm ratio")
        fig.savefig(
            output_prefix + f"species_bias.original_vs_corrected.scatterplot{label}.{suffix}.svg"
            .replace("..", "."),
            dpi=300, bbox_inches="tight")

    return metrics_r2


def plot_metrics_lineplot(metrics, keys=['read', 'umi', 'gene'], tail=None, suffix="", by_group=None):
    print(f"# {time.asctime()} - Plotting metrics per cell.")
    if tail:
        metrics = metrics.tail(tail)
    n = len(keys)

    if by_group is not None:
        groups = metrics[by_group].dropna().unique()
        suffix += f"by_{by_group}"
    else:
        by_group = "dummy"
        groups = ["dummy"]
        metrics.loc[:, by_group] = groups[0]

    fig, axis = plt.subplots(n, 2, figsize=(2 * 6, n * 3), tight_layout=True)
    for i, metric in enumerate(keys):
        for group in groups:
            print(metric, group)
            d = metrics.loc[metrics[by_group] == group, metric].sort_values()
            rank = d.rank(ascending=False, method="average")
            axis[i, 0].plot(rank, d, rasterized=True, label=group if group != "dummy" else None)
            axis[i, 1].plot(rank / rank.shape[0], d / d.shape[0], rasterized=True, label=group if group != "dummy" else None)
            axis[i, 0].axvline(expected_cell_number, linestyle="--", color="grey")
        for ax in axis[i, :]:
            ax.loglog()
            ax.set_xlabel("Barcodes")
            ax.set_ylabel(metric.capitalize() + "s")
            if by_group != "dummy":
                ax.legend()
    fig.savefig(
        output_prefix + f"metrics_per_cell.lineplot.{suffix}.svg"
        .replace("..", "."),
        dpi=300, bbox_inches="tight")


def plot_metrics_distplot(metrics, keys=['read', 'umi', 'gene'], tail=None, suffix="", by_group=None):
    print(f"# {time.asctime()} - Plotting metrics per cell.")
    if tail:
        metrics = metrics.tail(tail)
    n = len(keys)

    if by_group is not None:
        groups = metrics[by_group].dropna().unique()
        suffix += f"by_{by_group}"
    else:
        by_group = "dummy"
        groups = ["dummy"]
        metrics.loc[:, by_group] = groups[0]

    fig, axis = plt.subplots(
        1, n, figsize=(n * 3, 1 * 3),
        tight_layout=True, sharex="col", squeeze=True)
    for i, metric in enumerate(keys):
        axis[i].set_yscale("log")
        for j, group in enumerate(groups):
            d = metrics.loc[metrics[by_group] == group, :]
            sns.distplot(
                np.log10(d[metric]), kde=False, ax=axis[i],
                label=group if group != "dummy" else None)
        axis[i].set_xlabel(metric.capitalize() + "s (log)")
        axis[i].set_ylabel("Barcodes")
        if by_group != "dummy":
            axis[i].legend()
    fig.savefig(
        output_prefix + f"metrics_per_cell.distplot.{suffix}.svg"
        .replace("..", "."),
        dpi=300, bbox_inches="tight")


def plot_efficiency(metrics, keys=["umi", "gene"], tail=25000, suffix="", by_group=None, colour_by=None):
    """
    `metrics` must contain a column "read".
    """
    print(f"# {time.asctime()} - Plotting efficiency metrics per cell.")
    if tail:
        metrics = metrics.tail(tail)

    if colour_by is None:
        colour_by = "efficiency"
        if colour_by not in metrics.columns:
            metrics.loc[:, colour_by] = metrics['umi'] / metrics['read']

    if by_group is not None:
        groups = metrics[by_group].dropna().unique()
        n = len(groups)
        suffix += f"by_{by_group}"
    else:
        by_group = "dummy"
        groups = ["dummy"]
        metrics.loc[:, by_group] = groups[0]
        n = 1

    fig, axis = plt.subplots(n, 2, figsize=(2 * 4, n * 4), tight_layout=True, squeeze=False)
    fig.suptitle(f"Experiment performance: {sample_name}", ha="center")
    for i, group in enumerate(groups):
        d = metrics.loc[metrics[by_group] == group, :]
        if group != "dummy":
            axis[i, 0].set_title(f"{by_group.capitalize()}: {group}")
        for j, metric in enumerate(keys):
            # # fit curve
            (m, b), pcov = scipy.optimize.curve_fit(lin_func, d["read"], d[metric])
            # lm = LinearRegression().fit(d['read'].values.reshape(-1, 1), d[metric])
            # assert np.allclose(m, lm.coef_)
            # assert np.allclose(b, lm.intercept_)
            axis[i, j].text(
                0.5e5,
                expected_cell_number,
                s="y = {:.3f}x + {:.1f}".format(m, b),
                ha="center")

            col = axis[i, j].scatter(
                d['read'],
                d[metric],
                c=d[colour_by], rasterized=True, alpha=0.2, s=2)
            add_colorbar_to_axis(col, label=colour_by)
            axis[i, j].loglog()
            axis[i, j].set_xlabel("Reads per cell")
            axis[i, j].set_ylabel(f"Useful {metric}s per cell")

            # Plot regression line
            x = np.linspace(d["read"].min(), d["read"].max(), num=1000)
            axis[i, j].plot(x, lin_func(x, m, b), color="orange")
            # Plot 10k annotation
            y = lin_func(1e4, m, b)
            axis[i, j].text(
                1e4, y,
                s=f"{metric.capitalize()}s recovered\nwith 10.000\nreads:\n{y:.2f}",
                ha="left")
            # Plot X == Y
            xmax = d["read"].max()
            xmax += xmax * 0.1
            ymax = d[metric].max()
            ymax += ymax * 0.1
            x = np.linspace(0, ymax, num=2)
            y = lin_func(x, 1, 0)
            axis[i, j].plot(x, y, linestyle="--", color="black", linewidth=0.5)

            # Plot lines at relevant metrics
            for h in [100, 250, 500, 1000, expected_cell_number]:
                axis[i, j].axhline(h, linestyle="--", color="black", linewidth=0.5)
            for v in [10000, 100000]:
                axis[i, j].axvline(v, linestyle="--", color="black", linewidth=0.5)
    fig.savefig(
        output_prefix + f"performance_per_cell.scatter.coloured_by_{colour_by}.{suffix}.svg"
        .replace("..", "."),
        dpi=300, bbox_inches="tight")


def plot_species_mixing(metrics, norm=False, tail=None, suffix="", cmaps=None, zoom_in_area=3000, axislims=10000):
    print(f"# {time.asctime()} - Plotting species mixtures.")
    if tail is None:
        tail = 2 * expected_cell_number
    if tail is not False:
        metrics = metrics.tail(tail)

    if norm:
        suffix += "_norm"

    if cmaps is None:
        cmaps = [get_custom_cmap()]  #, plt.get_cmap("coolwarm"), plt.get_cmap("Spectral_r")]
    for attr, label, kwargs in [
            ('sp_ratio_norm' if norm else 'sp_ratio', 'coloured_by_ratio', {"vmin": 0, "vmax": 1}),
            # ('doublet_norm' if norm else 'doublet', 'coloured_by_doublet', {"vmin": -1, "vmax": 1}),
    ]:
        for cmap in cmaps:
            n_panels = 4
            fig, axis = plt.subplots(1, n_panels, figsize=(n_panels * 4, 4), tight_layout=True)
            for ax in axis[:3]:
                col = ax.scatter(
                    metrics['mouse'], metrics['human_norm' if norm else "human"],
                    c=metrics[attr], cmap=cmap,
                    s=2, alpha=0.25,
                    rasterized=True, **kwargs)
                ax.set_xlabel("Mouse (UMIs)")
                ax.set_ylabel("Human (UMIs)")
                add_colorbar_to_axis(col, label="Species fraction")
            v = metrics['total_norm' if norm else 'total'].max()
            v = min(axislims, v)
            axis[1].set_xlim((-(v / 30.), v))
            axis[1].set_ylim((-(v / 30.), v))
            axis[2].set_xlim((-(zoom_in_area / 30.), zoom_in_area))
            axis[2].set_ylim((-(zoom_in_area / 30.), zoom_in_area))
            axis[3].scatter(
                metrics['mouse'], metrics['human' if norm else "human"],
                c=metrics[attr], cmap=cmap,
                s=1, alpha=0.1,
                rasterized=True, **kwargs)
            axis[3].loglog()
            axis[3].set_xlabel("Mouse (UMIs)")
            axis[3].set_ylabel("Human (UMIs)")
            fig.savefig(
                output_prefix + f"species_mix.{suffix}.{label}.{cmap.name}.svg"
                .replace("..", "."),
                dpi=300, bbox_inches="tight")


def gather_stats_per_well(metrics, seq_content=False, save_intermediate=True):
    print(f"# {time.asctime()} - Gathering metrics per combinatorial indexing well.")
    # # # see how many of the round1 barcodes (well) are captured
    umis_per_well = metrics.groupby('r1')['umi'].sum().sort_values()
    umis_per_well.name = "umis"
    if save_intermediate:
        to_pickle(umis_per_well, 'umis_per_well')

    # # # see how many droplets per well
    droplets_per_well = metrics.reset_index().groupby('r1')['r2'].nunique().sort_values()
    droplets_per_well.name = "droplets"
    if save_intermediate:
        to_pickle(droplets_per_well, 'droplets_per_well')

    well_metrics = umis_per_well.to_frame().join(droplets_per_well)
    if save_intermediate:
        to_pickle(well_metrics, 'well_metrics')

    if seq_content:
        s = well_metrics.index.to_series().str
        well_metrics = well_metrics.assign(
            at_content=s.count("AT|TA") / 11,
            gc_content=s.count("GC|CG") / 11,
            a_content=s.count("A|T") / 11,
            c_content=s.count("C|G") / 11)
        if save_intermediate:
            to_pickle(well_metrics, 'well_metrics')

    return well_metrics


def get_exact_matches(metrics, r1_whitelist, r2_whitelist, save_intermediate=True, plot=True):
    print(f"# {time.asctime()} - Filtering for exact barcode matches.")

    # Note to self, if metrics has r1_annotations, one could simply do isnull() to get barcodes matching r1

    r1_match = metrics.index.get_level_values("r1").isin(r1_whitelist)
    r2_match = metrics.index.get_level_values("r2").isin(r2_whitelist)

    if save_intermediate:
        to_pickle(r1_match, "r1_match", array=False)
        to_pickle(r2_match, "r2_match", array=False)

    if plot:
        plot_barcode_match_fraction(r1_match, r2_match)
    metrics_filtered = metrics.loc[r1_match & r2_match]

    if save_intermediate:
        to_pickle(metrics_filtered, "metrics_filtered")
    return metrics_filtered


def get_exact_matches_droplet(metrics, r2_whitelist, save_intermediate=True, plot=True, suffix="_r2_only"):
    print(f"# {time.asctime()} - Filtering for exact barcode matches.")

    # Note to self, if metrics has r1_annotations, one could simply do isnull() to get barcodes matching r1
    r2_match = metrics.index.get_level_values("r2").isin(r2_whitelist)

    if save_intermediate:
        to_pickle(r2_match, "r2_match" + suffix, array=False)

    metrics_filtered = metrics.loc[r2_match]

    if save_intermediate:
        to_pickle(metrics_filtered, "metrics_filtered" + suffix)
    return metrics_filtered


def get_stats_per_droplet(metrics, doublet_threshold=0.85, save_intermediate=True):
    metrics_droplet = metrics.groupby("r2")['read', 'umi', 'gene', 'human', 'mouse', 'total', 'max'].sum()

    metrics_droplet = (
        metrics_droplet.assign(
            ratio=metrics_droplet["max"] / metrics_droplet['total'],
            sp_ratio=metrics_droplet["human"] / metrics_droplet['total'])
        .sort_values("total"))
    metrics_droplet = metrics_droplet.assign(doublet=(metrics_droplet['ratio'] < doublet_threshold).astype(int).replace(0, -1))

    # Assess species bias
    r = dict()
    for f in [0.1, 0.2, 0.25, 0.5, 0.75] + list(range(1, 1000, 10)) + [1.5]:
        t = metrics_droplet.tail(int(expected_cell_number * f))
        r[expected_cell_number * f] = t['mouse'].mean() / t['human'].mean()
    r = pd.Series(r).sort_index()

    # Add normalized metrics_droplet to stats
    metrics_droplet.loc[:, 'human_norm'] = metrics_droplet['human'] * r[int(expected_cell_number)]
    metrics_droplet.loc[:, 'total_norm'] = metrics_droplet[['mouse', 'human_norm']].sum(1)
    metrics_droplet.loc[:, 'max_norm'] = metrics_droplet[['mouse', 'human_norm']].max(1)
    metrics_droplet.loc[:, 'ratio_norm'] = metrics_droplet['max_norm'] / metrics_droplet['total_norm']
    metrics_droplet.loc[:, 'sp_ratio_norm'] = metrics_droplet['human_norm'] / metrics_droplet['total_norm']
    metrics_droplet.loc[:, 'doublet_norm'] = (metrics_droplet.loc[:, 'ratio_norm'] < doublet_threshold).astype(int).replace(0, -1)

    if save_intermediate:
        to_pickle(metrics_droplet, "metrics_droplet")

    return metrics_droplet


def plot_well_stats(well_metrics, keys=["droplets", "umis"], tail=None, suffix=""):
    print(f"# {time.asctime()} - Plotting metrics per combinatorial indexing well.")
    if tail is None:
        tail = n_round1_barcodes * 2
    well_metrics = well_metrics.tail(tail)

    n = len(keys)
    v = well_metrics[keys].min().min()
    fig, axis = plt.subplots(n, 1, figsize=(2 * 3, n * 3), tight_layout=True)
    for i, metric in enumerate(keys):
        d = well_metrics[metric].sort_values()
        rank = d.rank(ascending=False, method="average")
        axis[i].plot(rank, d, rasterized=True)
        axis[i].set_ylim((v - (v * 0.5), axis[i].get_ylim()[1] * 1.5))
        axis[i].set_yscale("log")
        axis[i].set_xlabel("Well")
        axis[i].set_ylabel(metric.capitalize())
    fig.savefig(
        output_prefix + f"metrics_per_well.lineplot.{suffix}.svg"
        .replace("..", "."),
        dpi=300, bbox_inches="tight")

    # # # # observe whether barcode sequence content relates with abundance
    if well_metrics.columns.str.endswith("_content").sum() < 1:
        return

    cols = well_metrics.columns[well_metrics.columns.str.endswith("_content")]
    fig, axis = plt.subplots(n, len(cols), figsize=(len(cols) * 3, n * 3), tight_layout=True)
    for j, col in enumerate(cols):
        for i, key in enumerate(keys):
            axis[i, j].scatter(well_metrics[col], well_metrics[key], rasterized=True, s=2, alpha=0.5)
            axis[i, j].set_xlim((-0.1, 1.1))
            axis[i, j].set_xlabel(col.replace("_", " ").upper())
            axis[i, j].set_ylabel(key.capitalize() + " per well")
    fig.savefig(
        output_prefix + f"metrics_per_well.sequence_content.{suffix}.svg"
        .replace("..", "."),
        dpi=300, bbox_inches="tight")


def plot_barcode_match_fraction(r1, r2, suffix=""):
    print(f"# {time.asctime()} - Plotting fraction of correct barcodes.")

    d = pd.DataFrame(
        [[r1.sum(), r2.sum(), (r1 & r2).sum()], [r1.shape[0]] * 3],
        index=['correct', 'total'], columns=['r1', 'r2', 'both'])
    d.loc["fraction of exact match", :] = d.loc['correct'] / d.loc['total']

    grid = sns.catplot(
        data=d.reset_index().melt(id_vars="index"),
        x="variable", y="value", col="index",
        sharey=False, kind="bar", height=3, aspect=1)
    grid.fig.savefig(
        output_prefix + f"barcode_matching.barplot.{suffix}.svg"
        .replace("..", "."),
        dpi=300, bbox_inches="tight")


def plot_comparison_to_10x(c, d, suffix=""):
    # # Compare side by side in logster
    for attr, label, kwargs in [
            ('sp_ratio', 'coloured_by_ratio', {"cmap": get_custom_cmap(), "vmin": 0, "vmax": 1}),
            # ('doublet', 'coloured_by_doublet', {"cmap": "tab20b", "vmin": -1, "vmax": 1}),
            # ('total', 'coloured_by_total', {"vmin": 100, "vmax": 10000}),
    ]:
        fig, axis = plt.subplots(2, 3, figsize=(3 * 4, 2 * 4), tight_layout=True)
        for ax in axis[:, 0]:
            ax.set_title("10X")
            ax.scatter(
                d['mouse'], d['human'],
                c=d[attr],
                s=1, alpha=0.1,
                rasterized=True, **kwargs)
        for ax in axis[:, 1]:
            ax.set_title("scifi-RNA-seq")
            ax.scatter(
                c['mouse'], c['human'],
                c=c[attr],
                s=1, alpha=0.1,
                rasterized=True, **kwargs)
        for ax in axis[:, 2]:
            ax.set_title("scifi-RNA-seq coloured as 10X only")
            c.loc[:, "10X_attr"] = d.loc[c.index.get_level_values(1), attr].values
            ax.scatter(
                c['mouse'], c['human'],
                c=c["10X_attr"],
                s=1, alpha=0.05,
                rasterized=True, **kwargs)
        for ax in axis.flatten():
            ax.set_xlabel("Mouse (UMIs)")
            ax.set_ylabel("Human (UMIs)")

        for ax in axis[0, :]:
            ax.set_xlim((-(7500 / 30.), 7500))
            ax.set_ylim((-(7500 / 30.), 7500))
        for ax in axis[1, :]:
            ax.loglog()
        fig.savefig(
            output_prefix + f"species_mix.comparison_to_droplet.{label}.{suffix}.svg"
            .replace("..", "."),
            dpi=300, bbox_inches="tight")

    fig, axis = plt.subplots(3, 3, figsize=(3 * 3, 3 * 3), tight_layout=True)
    attr, label, kwargs = ('doublet', 'coloured_by_doublet', {"cmap": "tab20b", "vmin": 0, "vmax": 1})
    for i, t in enumerate(pd.np.arange(0.6, 0.8, 0.1)):
        # get only duplets from 10X
        p_10x = d.loc[d.loc[:, 'ratio'] < t, :]
        # get scifi cells which were in those duplet droplets
        p_spc = c.loc[c.index.get_level_values(1).isin(p_10x.index)]
        axis[i, 0].set_title(f"10X\n(duplets only at {t:.2f} purity)")
        axis[i, 0].scatter(
            p_10x['mouse'], p_10x['human'],
            c=p_10x[attr] + 1,
            s=1, alpha=0.1,
            rasterized=True, **kwargs)
        axis[i, 1].set_title(f"scifi-RNA-seq\n(only droplets which are duplets in 10X at {t:.2f} purity)")
        axis[i, 1].scatter(
            p_spc['mouse'], p_spc['human'],
            c=p_spc[attr] + 1,
            s=1, alpha=0.1,
            rasterized=True, **kwargs)
        axis[i, 2].set_title(f"scifi-RNA-seq\n(only droplets which are duplets in 10X at {t:.2f} purity)")
        axis[i, 2].scatter(
            p_10x['mouse'], p_10x['human'],
            color=plt.get_cmap('Pastel1')(0),
            s=1, alpha=0.1,
            rasterized=True, label="10X")
        axis[i, 2].scatter(
            p_spc['mouse'], p_spc['human'],
            color=plt.get_cmap('Pastel1')(1),
            s=1, alpha=0.1,
            rasterized=True, label="scifi-RNA-seq")
        axis[i, 2].legend()

    for ax in axis.flatten():
        ax.set_xlim((-(3000 / 30.), 3000))
        ax.set_ylim((-(3000 / 30.), 3000))

    fig.savefig(
        output_prefix + f"species_mix.comparison_to_only_10X.only_doublets.{label}.{suffix}.svg"
        .replace("..", "."),
        dpi=300, bbox_inches="tight")


def to_pickle(obj, name, array=True, only_array=False):
    print(f"# {time.asctime()} - Saving {name} to pickle.")
    if array:
        pickle.dump(obj.values, open(input_prefix + f"{name}.values.pickle", 'wb'), protocol=-1)
    if only_array:
        return
    pickle.dump(obj, open(input_prefix + f"{name}.pickle", 'wb'), protocol=-1)


def from_pickle(key, array=False):
    print(f"# {time.asctime()} - Loading {key} from pickle.")
    if array:
        return pickle.load(open(input_prefix + f"{key}.values.pickle", 'rb'))
    return pickle.load(open(input_prefix + f"{key}.pickle", 'rb'))

    # df = pickle.load(open(input_prefix + "all.pickle", 'rb'))


def cells_per_droplet_stats(cells_per_droplet):
    # Observe statistical properties
    def poisson(k, lamb):
        from scipy.special import factorial
        return np.exp(-lamb) * (lamb ** k / factorial(k))

    def log_linear_poisson(k, lamb):
        from scipy.special import factorial
        return np.log(np.exp(-lamb) * (lamb ** k / factorial(k)))

    cells_per_droplet_counts = cells_per_droplet.value_counts().sort_index()

    lamb, cov_matrix = scipy.optimize.curve_fit(
        log_linear_poisson,
        cells_per_droplet_counts.index.astype(int),
        cells_per_droplet_counts)

    print(lamb)

    fig, axis = plt.subplots(1, 2, figsize=(3 * 2, 3), tight_layout=True)
    bins = cells_per_droplet_counts.shape[0]
    for ax in axis:
        sns.distplot(cells_per_droplet, bins=bins, kde=False, ax=ax, label="Real")
    x = np.arange(0, cells_per_droplet_counts.index.max())
    y_hat = scipy.stats.poisson(lamb).pmf(x)
    y_hat *= cells_per_droplet.shape[0]
    for ax in axis:
        ax.plot(x + 0.5, y_hat)  # the 0.5 is just to center on the middle of the histogram bins

    cpd = scipy.stats.poisson(lamb).rvs(cells_per_droplet.shape[0])
    for ax in axis:
        sns.distplot(cpd, bins=bins, kde=False, ax=ax, norm_hist=False, label="Poisson")
    for ax in axis:
        ax.set_xlabel("Cells")
        ax.set_ylabel("Droplets")
        ax.legend()
    axis[1].set_yscale("log")
    axis[1].set_ylim(bottom=0.1)
    fig.savefig(
        output_prefix + f"cells_per_droplet.poissonian_properties.svg",
        dpi=300, bbox_inches="tight")


def add_colorbar_to_axis(collection, label=None, position="right", size="5%", pad=0.05):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(collection.axes)
    cax = divider.append_axes(position, size=size, pad=pad)
    cbar = plt.colorbar(mappable=collection, cax=cax, label=label, alpha=1)
    cbar.solids.set_edgecolor("face")
    cbar.solids.set_rasterized(True)


def lin_func(x, m, b):
    return m * x + b


def get_custom_cmap(vmin=-1.0, vmax=1.0, inner=0.4):
    norm = matplotlib.colors.Normalize(-1, 1)
    colors = [
        [norm(vmin), "cornflowerblue"],
        [norm(-inner), "crimson"],
        [norm(inner), "crimson"],
        [norm(vmax), "burlywood"]]
    cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
    cmap.name = "custom"
    return cmap


def parallel_groupby_apply(df, levels, function):
    from joblib import Parallel, delayed
    import multiprocessing

    g = df.groupby(levels)
    res = Parallel(n_jobs=multiprocessing.cpu_count())(delayed(function)(group) for name, group in g)
    return pd.concat(res)


if __name__ == "__main__":
    sns.set(context="paper", style="ticks", palette="colorblind", color_codes=True)
    matplotlib.rcParams["svg.fonttype"] = "none"
    # Don't use LaTeX for rendering
    matplotlib.rcParams["text.usetex"] = False
    sys.exit(main())
