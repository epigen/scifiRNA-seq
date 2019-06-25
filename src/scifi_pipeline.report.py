#! /bin/env python

import sys
import os
from argparse import ArgumentParser
import time
import pickle

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import scipy


def parse_args():
    parser = ArgumentParser()
    parser.add_argument(
        "--metric-file", dest="metric_files",
        help="Input metric files from the output of sci-RNA.summarizer.py. Only one can be provided.",
        required=True, nargs="+")
    parser.add_argument(
        "--output-prefix", dest="output_prefix",
        help="Absolute path prefix for output files. ",
        required=True)
    parser.add_argument(
        "--plotting-attributes", dest="plotting_attributes", nargs="*",
        help="Additional attributes to plot metrics grouped by. "
             "Several can be given. Example: 'donor_id sex plate'.")
    parser.add_argument(
        "--plate-column", dest="plate_column", default="plate",
        help="Name of column in input metric files containing r1 plate information. "
             "Defaults to 'plate'.")
    parser.add_argument(
        "--well-column", dest="well_column", default="plate_well",
        help="Name of column in input metric files containing r1 well plate information. "
             "Defaults to 'plate_well'.")
    parser.add_argument(
        "--droplet-column", dest="droplet_column", default="r2",
        help="Name of column in input metric files containing r2 (droplet) information. "
             "Defaults to 'r2'.")
    parser.add_argument(
        "--nrows", dest="nrows",
        help="An optional number of rows of input to process. "
             "Useful for testing.")
    parser.add_argument(
        "--species-mixture", dest="species_mixture", action="store_true",
        help="Whether the experiment is a species mixture. "
             "Default is false.")
    parser.add_argument(
        "--expected-cell-number", dest="expected_cell_number", default=200000, type=int,
        help="Expected experiment yield number.")
    parser.add_argument(
        "--save-intermediate", dest="save_intermediate", action="store_true",
        help="Whether to save intermediate pickle files.")

    # # Example:
    # args = parser.parse_args([
        # "--metric-file", "/scratch/lab_bock/shared/projects/sci-rna/data/SCI024_Tcell_s/SCI024_Tcell_s.exon.metrics.csv.gz",
        # "--output-prefix", "results/SCI024_Tcell_s.exon.",
        # "--plotting-attributes", "plate", "donor_id"])
    args = parser.parse_args()


    if not args.output_prefix.endswith("."):
        args.output_prefix += "."

    # expand glob if needed
    if len(args.metric_files) == 1:
        if "*" in args.metric_files[0]:
            args.metric_files = glob(args.metric_files[0])

    if not os.path.exists(os.path.dirname(args.output_prefix)):
        os.makedirs(os.path.dirname(args.output_prefix))

    return args


def main():
    global args
    args = parse_args()

    # read text files
    print(args.metric_files)
    metrics = load_metrics(args.metric_files, nrows=args.nrows)

    # # Plot
    # # # Loglog line plot of rank vs abundance
    t = args.expected_cell_number * 5  # number of barcodes to plot
    plot_metrics_lineplot(metrics, tail=t)  # knee-plot
    for attribute in args.plotting_attributes:
        plot_metrics_lineplot(metrics, tail=t, by_group=attribute)  # e.g. used for fixation, no fixation

    # # # Efficiency plot (UMI and gene yield in relation to reads per cell)
    t = int(args.expected_cell_number * 5)
    plot_efficiency(metrics, tail=t)
    plot_efficiency(metrics, tail=t, colour_by="unique_fraction")

    # Well and droplet inspection
    well_metrics = gather_stats_per_well(metrics, seq_content=True)
    plot_well_stats(well_metrics, tail=None, suffix="")

    cells_per_droplet = metrics.groupby(args.droplet_column)[args.well_column].nunique()
    cells_per_droplet_stats(cells_per_droplet)

    # Write expression matrix in sparse format as H5AD
    if len(args.metric_files) == 1:
        print(f"# {time.asctime()} - Reading sparse expression matrix.")
        expr = pd.read_csv(args.metric_files[0].replace("metrics", "expression"))
        write_gene_expression_matrix(expr)


def write_gene_expression_matrix(expr, output_file=None, file_format="h5ad"):
    from scipy.sparse import csr_matrix
    import anndata as an

    if file_format not in ['h5ad', 'loom']:
        raise ValueError("Output format must be one of 'h5ad' or 'loom'.")

    n = expr['umi'].isnull().sum()
    if n > 0:
        print(f"# {time.asctime()} - Dropping {n} entries with null UMI values.")
        expr = expr.dropna(subset=['umi'])

    print(f"# {time.asctime()} - Adding cell label.")
    expr = expr.assign(cell=expr[args.plate_column] + '-' + expr[args.well_column] + '-' + expr[args.droplet_column])
    print(f"# {time.asctime()} - Factorizing cells and genes.")
    cell, unique_cell = expr['cell'].factorize()
    gene, unique_gene = expr['gene'].factorize()

    print(f"# {time.asctime()} - Creating AnnData object.")
    a = an.AnnData(
        csr_matrix((expr['umi'], (cell, gene)), shape=(len(unique_cell), len(unique_gene))),
        dtype=np.int)
    a.obs.index = unique_cell
    a.var.index = unique_gene

    print(f"# {time.asctime()} - Writing h5ad object to disk.")
    if output_file is None:
        output_file = args.output_prefix + file_format
    if file_format == "h5ad":
        a.write(output_file)
    if file_format == "loom":
        a.write_loom(output_file)


def load_metrics(files, **kwargs):
    print(f"# {time.asctime()} - Reading files.")
    pieces = list()
    for file in files:
        print(f"# {time.asctime()} - Reading file {file}.")
        d = pd.read_csv(
            file,
            error_bad_lines=False, engine='c',
            compression='gzip', **kwargs)
        pieces.append(d)
        print(f"# {time.asctime()} - Done with file {file}. {d.shape[0]} lines.")

    print(f"# {time.asctime()} - Concatenating parts.")
    return pd.concat(pieces).sort_values("umi")


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
        t = metrics_r2.tail(int(args.expected_cell_number * f))
        r[args.expected_cell_number * f] = t['mouse'].mean() / t['human'].mean()
    r = pd.Series(r).sort_index()

    # Add normalized metrics to stats
    metrics_r2.loc[:, 'human_norm'] = metrics_r2['human'] * r[int(args.expected_cell_number)]
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
    axis.axvline(args.expected_cell_number, linestyle="--", color="grey")
    axis.set_xlabel("Top N barcodes taken into account")
    axis.set_ylabel("Ratio of mean UMIs per cell\nbetween mouse and human")
    axis.set_xscale("log")
    fig.savefig(
        args.output_prefix + f"species_bias.lineplot.{suffix}.svg"
        .replace("..", "."),
        dpi=300, bbox_inches="tight")

    # Plot illustration of normalization procedure
    metrics2 = metrics_r2.tail(int(10 * args.expected_cell_number))

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
            ax.plot((0, v), (0, v * (1 / r[int(args.expected_cell_number)])), linestyle="--", color="orange")
            ax.set_ylabel("Human (UMIs)")

        for ax in axis.flatten():
            ax.set_xlabel("Mouse (UMIs)")
        axis[0, 0].set_title("Coloured by ratio")
        axis[0, 1].set_title("Coloured by norm ratio")
        fig.savefig(
            args.output_prefix + f"species_bias.original_vs_corrected.scatterplot{label}.{suffix}.svg"
            .replace("..", "."),
            dpi=300, bbox_inches="tight")

    return metrics_r2


def plot_metrics_lineplot(metrics, keys=['read', 'umi', 'gene'], tail=None, suffix="", by_group=None, always_legend=False):
    def min_max(x):
        return (x - x.min()) / (x.max() - x.min())

    print(f"# {time.asctime()} - Plotting metrics per cell.")
    if tail:
        metrics = metrics.tail(tail)
    n = len(keys)

    if by_group is not None:
        groups = sorted(metrics[by_group].dropna().unique())
        suffix += f"by_{by_group}"
    else:
        by_group = "dummy"
        groups = ["dummy"]
        metrics.loc[:, by_group] = groups[0]

    colors = sns.color_palette('husl', n_colors=len(groups))
    styles = ['solid', 'dashed', 'dashdot', 'dotted']
    num_styles = len(styles)

    fig, axis = plt.subplots(n, 2, figsize=(2 * 6, n * 3))#, tight_layout=True)
    for i, metric in enumerate(keys):
        for j, group in enumerate(groups):
            d = metrics.loc[metrics[by_group] == group, metric].sort_values()
            d_norm = min_max(d)
            rank = d.rank(ascending=False, method="average", pct=False)
            rank_norm = min_max(rank)
            kwargs = {
                "color": colors[j], "linestyle": styles[j % num_styles],
                "rasterized": True, "label": group if group != "dummy" else None}
            axis[i, 0].plot(rank, d, **kwargs)
            axis[i, 1].plot(rank_norm, d_norm, **kwargs)
        axis[i, 0].axvline(args.expected_cell_number, linestyle="--", color="grey")
        for l, ax in enumerate(axis[i, :]):
            ax.loglog()
            ax.set_xlabel("Barcodes")
            ax.set_ylabel(metric.capitalize() + "s")
            # add legend only to one panel if requested
            if (by_group != "dummy") and (always_legend or ((i == 0) and (l == 1))):
                ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", borderaxespad=0., ncol=1)
                # ax.legend_.set_in_layout(False)
    fig.savefig(
        args.output_prefix + f"metrics_per_cell.lineplot.{suffix}.svg"
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

    colors = sns.color_palette('husl', n_colors=len(groups))

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
        args.output_prefix + f"metrics_per_cell.distplot.{suffix}.svg"
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

    fig, axis = plt.subplots(
        n, 2, figsize=(2 * 4, n * 4),
        tight_layout=True, squeeze=False, sharex="col", sharey="col")
    fig.suptitle(f"Experiment performance: {args.output_prefix}", ha="center")
    for i, group in enumerate(groups):
        d = metrics.loc[metrics[by_group] == group, :]
        if group != "dummy":
            axis[i, 0].set_title(f"{by_group.capitalize()}: {group}")
        for j, metric in enumerate(keys):
            # # fit curve
            try:
                (m, b), pcov = scipy.optimize.curve_fit(lin_func, d["read"], d[metric])
                # lm = LinearRegression().fit(d['read'].values.reshape(-1, 1), d[metric])
                # assert np.allclose(m, lm.coef_)
                # assert np.allclose(b, lm.intercept_)
                axis[i, j].text(
                    0.5e5,
                    args.expected_cell_number,
                    s="y = {:.3f}x + {:.1f}".format(m, b),
                    ha="center")
                # Plot regression line
                x = np.linspace(d["read"].min(), d["read"].max(), num=1000)
                axis[i, j].plot(x, lin_func(x, m, b), color="orange")
            except TypeError:
                print(f"Couldn't find a fit for {group}, {metric}.")
                pass

            # Plot cells
            col = axis[i, j].scatter(
                d['read'],
                d[metric],
                c=d[colour_by], rasterized=True, alpha=0.2, s=2)
            add_colorbar_to_axis(col, label=colour_by)
            axis[i, j].loglog()
            axis[i, j].set_xlabel("Reads per cell")
            axis[i, j].set_ylabel(f"Useful {metric}s per cell")

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
            for h in [100, 250, 500, 1000, args.expected_cell_number]:
                axis[i, j].axhline(h, linestyle="--", color="black", linewidth=0.5)
            for v in [10000, 100000]:
                axis[i, j].axvline(v, linestyle="--", color="black", linewidth=0.5)
    fig.savefig(
        args.output_prefix + f"performance_per_cell.scatter.coloured_by_{colour_by}.{suffix}.svg"
        .replace("..", "."),
        dpi=300, bbox_inches="tight")


def plot_species_mixing(metrics, norm=False, tail=None, suffix="", cmaps=None, zoom_in_area=3000, axislims=10000):
    print(f"# {time.asctime()} - Plotting species mixtures.")
    if tail is None:
        tail = 2 * args.expected_cell_number
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
                args.output_prefix + f"species_mix.{suffix}.{label}.{cmap.name}.svg"
                .replace("..", "."),
                dpi=300, bbox_inches="tight")


def gather_stats_per_well(metrics, seq_content=False, save_intermediate=True):
    print(f"# {time.asctime()} - Gathering metrics per combinatorial indexing well.")
    # # # see how many of the round1 barcodes (well) are captured
    umis_per_well = metrics.groupby(args.well_column)['umi'].sum().sort_values()
    umis_per_well.name = "umis"
    if save_intermediate:
        to_pickle(umis_per_well, 'umis_per_well')

    # # # see how many droplets per well
    droplets_per_well = metrics.reset_index().groupby(args.well_column)['r2'].nunique().sort_values()
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
        plot_umi_match_fraction(metrics['umi'], r1_match, r2_match)
    metrics_filtered = metrics.loc[r1_match & r2_match]

    if save_intermediate:
        to_pickle(metrics_filtered, "metrics_filtered")
    return metrics_filtered

    # from Levenshtein import distance
    # from joblib import Parallel, delayed

    # query = metrics.index.get_level_values("r1")[~r1_match][:100]
    # reference = r1_whitelist
    # distances = Parallel(n_jobs=-1)(
    #     delayed(distance)(q, r)
    #     for q in query
    #     for r in reference)

    metrics = metrics.assign(r1_match=r1_match)
    m = metrics.tail(2000000).groupby('r1_match').mean()
    mm = m.reset_index().melt(id_vars=['r1_match'])
    grid = sns.catplot(
        data=mm, x="r1_match", col="variable", y="value", kind='bar',
        sharey=False, col_wrap=4, height=2, aspect=1)
    grid.savefig(args.output_prefix + "metrics_dependent_on_r1_match.svg", bbox_inches="tight", dpi=300)


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
        t = metrics_droplet.tail(int(args.expected_cell_number * f))
        r[args.expected_cell_number * f] = t['mouse'].mean() / t['human'].mean()
    r = pd.Series(r).sort_index()

    # Add normalized metrics_droplet to stats
    metrics_droplet.loc[:, 'human_norm'] = metrics_droplet['human'] * r[int(args.expected_cell_number)]
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
    if tail is not None:
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
        args.output_prefix + f"metrics_per_well.lineplot.{suffix}.svg"
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
        args.output_prefix + f"metrics_per_well.sequence_content.{suffix}.svg"
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
        args.output_prefix + f"barcode_matching.barplot.{suffix}.svg"
        .replace("..", "."),
        dpi=300, bbox_inches="tight")


def plot_umi_match_fraction(umi, r1_match, r2_match, suffix=""):
    print(f"# {time.asctime()} - Plotting fraction of correct barcodes.")

    d = pd.DataFrame(
        [[umi.loc[r1_match].sum(), umi.sum()], [umi.loc[r2_match].sum(), umi.sum()]],
        index=['correct', 'total'], columns=['r1', 'r2'])
    d.loc["fraction of exact match", :] = d.loc['correct'] / d.loc['total']

    grid = sns.catplot(
        data=d.reset_index().melt(id_vars="index"),
        x="variable", y="value", col="index",
        sharey=False, kind="bar", height=3, aspect=1)
    grid.fig.savefig(
        args.output_prefix + f"umi_matching.barplot.{suffix}.svg"
        .replace("..", "."),
        dpi=300, bbox_inches="tight")


def sequence_content_analysis(metrics, whitelist, barcode="r1", head=5000000):
    from Bio import motifs
    from Bio.Seq import Seq

    r1 = metrics.head(head).index.get_level_values(barcode).to_series()

    in_ = r1.isin(whitelist.tolist())
    r1_right = r1[in_]
    r1_wrong = r1[~in_]

    motif_right = motifs.create([Seq(x) for x in r1_right if "N" not in x])
    motif_wrong = motifs.create([Seq(x) for x in r1_wrong if "N" not in x])

    r = pd.DataFrame(motif_right.pwm)
    w = pd.DataFrame(motif_wrong.pwm)
    c = np.log2(w / r)

    fig, axis = plt.subplots(1, 3, figsize=(3 * 3, 3))
    kwargs = {"square": True}
    sns.heatmap(r.T, ax=axis[0], **kwargs)
    axis[0].set_title("Correct")
    sns.heatmap(w.T, ax=axis[1], **kwargs)
    axis[1].set_title("Wrong")
    kwargs = {"square": True, "cmap": "RdBu_r", "center": 0, "robust": True}
    sns.heatmap(c.T, ax=axis[2], **kwargs)
    axis[2].set_title("Difference")
    fig.savefig(args.output_prefix + f"barcode_{barcode}.sequence_content.svg", bbox_inches="tight", dpi=300)


def to_pickle(obj, name, array=True, only_array=False):
    print(f"# {time.asctime()} - Saving {name} to pickle.")
    if array:
        pickle.dump(obj.values, open(args.output_prefix + f"{name}.values.pickle", 'wb'), protocol=-1)
    if only_array:
        return
    pickle.dump(obj, open(args.output_prefix + f"{name}.pickle", 'wb'), protocol=-1)


def from_pickle(key, array=False):
    print(f"# {time.asctime()} - Loading {key} from pickle.")
    if array:
        return pickle.load(open(args.output_prefix + f"{key}.values.pickle", 'rb'))
    return pickle.load(open(args.output_prefix + f"{key}.pickle", 'rb'))

    # df = pickle.load(open(args.output_prefix + "all.pickle", 'rb'))


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

    fig, axis = plt.subplots(1, 2, figsize=(3 * 2, 3), tight_layout=True)
    # bins = cells_per_droplet_counts.shape[0]
    bins = None
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
        args.output_prefix + f"cells_per_droplet.poissonian_properties.svg",
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
