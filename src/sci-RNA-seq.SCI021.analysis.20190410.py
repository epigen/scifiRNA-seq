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
    annotation_file = os.path.join("metadata", "sciRNA-seq.SCI017.oligos_2019-02-11.csv")
    annotation = pd.read_csv(annotation_file)

    droplet_barcode_file = os.path.join("metadata", "737K-cratac-v1.reverse_complement.csv")
    droplets = pd.read_csv(droplet_barcode_file)

    # sample annotations
    # sample_name = "sci-RNA-seq_SCI020_3uL_reseq_4K"
    # n_parts = 4
    # n_round1_barcodes = 96
    # round2_barcodes = "reverse_complement"
    # expected_cell_number = 4000

    sample_name = "sci-RNA-seq_SCI021_125K"
    n_parts = 2
    n_round1_barcodes = 96
    round2_barcodes = "original"
    expected_cell_number = 125000

    # convenience
    sample_dir = os.path.join("data", sample_name)
    results_dir = "results"
    input_prefix = os.path.join(sample_dir, f"{sample_name}.")
    output_prefix = os.path.join(results_dir, f"{sample_name}.")

    # read text files
    df = load_data(sample_name, n_parts)
    print(f"# {time.asctime()} - Writing data to pickle.")
    pickle.dump(df, open(input_prefix + "all.pickle", 'wb'), protocol=-1)
    # df = pickle.load(open(input_prefix + ".all.pickle", 'rb'))

    # Gather metrics per cell
    metrics = gather_stats_per_cell(df)

    # # Plot
    # # # Loglog line plot of rank vs abundance
    plot_metrics_lineplot(metrics.tail(int(expected_cell_number * 10)))
    # # # Efficiency plot (UMI and gene yield in relation to reads per cell)
    plot_efficiency(metrics.tail(int(expected_cell_number * 2)), tail=False)
    # # # Species mixing - only for twice the number of expected cells
    plot_species_mixing(metrics.tail(int(expected_cell_number * 1.5)))

    # Well and droplet inspection
    well_metrics = gather_stats_per_well(metrics, seq_content=True)
    plot_well_stats(well_metrics, tail=None, suffix="")

    #

    # Now do the same only for exactly matching barcodes
    r1_ = metrics.index.get_level_values("r1").isin(annotation['barcode_sequence'])
    r2_ = metrics.index.get_level_values("r2").isin(droplets[round2_barcodes])
    plot_barcode_match_fraction(r1_, r2_)
    metrics_fil = metrics.loc[r1_ & r2_]

    plot_metrics_lineplot(metrics_fil.tail(expected_cell_number * 10), suffix="exact_match")
    # # # Efficiency plot (UMI and gene yield in relation to reads per cell)
    plot_efficiency(metrics_fil.tail(expected_cell_number * 5), tail=False, suffix="exact_match")
    # # # Species mixing - only for twice the number of expected cells
    plot_species_mixing(metrics_fil.tail(int(expected_cell_number * 1.5)), tail=False, suffix="exact_match")

    # Well and droplet inspection
    well_metrics_fil = gather_stats_per_well(metrics_fil, seq_content=True, save_intermediate=True)
    plot_well_stats(well_metrics_fil, tail=None, suffix="exact_match")


def gather_stats_per_cell(df, doublet_threshold=0.85, save_intermediate=True):
    print(f"# {time.asctime()} - Gathering metrics per cell.")

    # performance metrics
    reads_per_cell = df.groupby(['r1', 'r2'], sort=False)['read'].nunique().sort_values()
    if save_intermediate:
        to_pickle(reads_per_cell, 'reads_per_cell')
    umi_counts = df.groupby(['r1', 'r2', 'gene', 'pos'], sort=False)['umi'].nunique()
    if save_intermediate:
        to_pickle(umi_counts, 'umi_counts')
    umis_per_cell = umi_counts.groupby(['r1', 'r2'], sort=False).sum().sort_values()
    if save_intermediate:
        to_pickle(umis_per_cell, 'umis_per_cell')
    genes_per_cell = umi_counts.reset_index(level='gene').groupby(['r1', 'r2'], sort=False)['gene'].nunique().sort_values()
    if save_intermediate:
        to_pickle(genes_per_cell, 'genes_per_cell')

    # species mixing
    print(f"# {time.asctime()} - Gathering species-specific metrics per cell.")
    umi_counts2 = umi_counts.reset_index()
    umi_counts2 = (
        umi_counts2
        .assign(
            species=(
                umi_counts2['gene'].str.startswith("ENSG")
                .replace(True, "human").replace(False, "mouse"))))
    species_counts = umi_counts2.groupby(['r1', 'r2', 'species'])['umi'].sum()

    spc = species_counts.reset_index().pivot_table(
        index=['r1', 'r2'], columns='species', values='umi', fill_value=0)
    spc += 1
    spc = spc.assign(total=spc.sum(1), max=spc.max(1))
    spc = (
        spc.assign(
            ratio=spc["max"] / spc['total'],
            sp_ratio=spc["human"] / spc['total'])
        .sort_values("total"))
    spc = spc.assign(doublet=(spc['ratio'] < doublet_threshold).astype(int).replace(0, -1))
    if save_intermediate:
        to_pickle(spc, "spc")

    print(f"# {time.asctime()} - Joining all metrics.")
    metrics = reads_per_cell.to_frame().join(umis_per_cell).join(genes_per_cell).join(spc)
    if save_intermediate:
        to_pickle(metrics, "metrics")
    return metrics


def plot_metrics_lineplot(metrics, keys=['read', 'umi', 'gene'], tail=None, suffix=""):
    print(f"# {time.asctime()} - Plotting metrics per cell.")
    if tail:
        metrics = metrics.tail(tail)
    n = len(keys)
    fig, axis = plt.subplots(n, 1, figsize=(6, n * 3), tight_layout=True)
    for i, metric in enumerate(keys):
        d = metrics[metric].sort_values()
        rank = d.rank(ascending=False, method="average")
        axis[i].plot(rank, d, rasterized=True)
        axis[i].axvline(expected_cell_number, linestyle="--", color="grey")
        axis[i].loglog()
        axis[i].set_xlabel("Barcodes")
        axis[i].set_ylabel(metric.capitalize() + "s")
    fig.savefig(
        output_prefix + f"metrics_per_cell.lineplot.{suffix}.svg"
        .replace("..", "."),
        dpi=300, bbox_inches="tight")


def plot_efficiency(metrics, keys=["umi", "gene"], tail=25000, suffix=""):
    """
    `metrics` must contain a column "read".
    """
    print(f"# {time.asctime()} - Plotting efficiency metrics per cell.")
    if tail:
        metrics = metrics.tail(tail)

    if "efficiency" not in metrics.columns:
        metrics.loc[:, "efficiency"] = metrics['umi'] / metrics['read']

    ## TODO: plot per material
    fig, axis = plt.subplots(1, 2, figsize=(2 * 4, 1 * 4), tight_layout=True)
    fig.suptitle(f"Experiment efficiency: {sample_name}", ha="center")
    for i, metric in enumerate(keys):
        # # fit curve
        (m, b), pcov = scipy.optimize.curve_fit(f, metrics["read"], metrics[metric])
        # lm = LinearRegression().fit(metrics['read'].values.reshape(-1, 1), metrics[metric])
        # assert np.allclose(m, lm.coef_)
        # assert np.allclose(b, lm.intercept_)
        axis[i].text(
            0.5e5,
            expected_cell_number,
            s="y = {:.3f}x + {:.1f}".format(m, b),
            ha="center")

        col = axis[i].scatter(
            metrics['read'],
            metrics[metric],
            c=metrics['efficiency'], rasterized=True, alpha=0.2, s=2)
        add_colorbar_to_axis(col, label="Relative efficiency")
        axis[i].loglog()
        axis[i].set_xlabel("Reads per cell")
        axis[i].set_ylabel(f"Useful {metric}s per cell")

        # Plot regression line
        x = np.linspace(metrics["read"].min(), metrics["read"].max(), num=1000)
        axis[i].plot(x, f(x, m, b), color="orange")
        # Plot 10k annotation
        y = f(1e4, m, b)
        axis[i].text(
            1e4, y,
            s=f"{metric.capitalize()}s recovered\nwith 10.000\nreads:\n{y:.2f}",
            ha="left")
        # Plot X == Y
        xmax = metrics["read"].max()
        xmax += xmax * 0.1
        ymax = metrics[metric].max()
        ymax += ymax * 0.1
        x = np.linspace(0, ymax, num=2)
        y = f(x, 1, 0)
        axis[i].plot(x, y, linestyle="--", color="black", linewidth=0.5)

        # Plot lines at relevant metrics
        for h in [100, 250, 500, 1000, expected_cell_number]:
            axis[i].axhline(h, linestyle="--", color="black", linewidth=0.5)
        for v in [10000, 100000]:
            axis[i].axvline(v, linestyle="--", color="black", linewidth=0.5)
    fig.savefig(
        output_prefix + f"performance_per_cell.scatter.{suffix}.svg"
        .replace("..", "."),
        dpi=300, bbox_inches="tight")


def plot_species_mixing(metrics, tail=25000, suffix="", cmaps=None, zoom_in_area=3000, axislims=10000):
    print(f"# {time.asctime()} - Plotting species mixtures.")
    if tail:
        metrics = metrics.tail(tail)
    if cmaps is None:
        cmaps = [plt.get_cmap("Spectral_r"), plt.get_cmap("coolwarm"), get_custom_cmap()]
    for attr, label, kwargs in [
            ('doublet', 'coloured_by_doublet', {"vmin": -1, "vmax": 1}),
            ('sp_ratio', 'coloured_by_ratio', {"vmin": 0, "vmax": 1}),
    ]:
        for cmap in cmaps:
            n_panels = 4
            fig, axis = plt.subplots(1, n_panels, figsize=(n_panels * 4, 4), tight_layout=True)
            for ax in axis[:3]:
                col = ax.scatter(
                    metrics['mouse'], metrics['human'],
                    c=metrics[attr], cmap=cmap,
                    s=2, alpha=0.25,
                    rasterized=True, **kwargs)
                ax.set_xlabel("Mouse (UMIs)")
                ax.set_ylabel("Human (UMIs)")
                add_colorbar_to_axis(col, label="Species fraction")
            v = metrics['total'].max()
            v = min(axislims, v)
            axis[1].set_xlim((-(v / 30.), v))
            axis[1].set_ylim((-(v / 30.), v))
            axis[2].set_xlim((-(zoom_in_area / 30.), zoom_in_area))
            axis[2].set_ylim((-(zoom_in_area / 30.), zoom_in_area))
            axis[3].scatter(
                metrics['mouse'], metrics['human'],
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


def load_data(sample_name, n_parts=2):
    print(f"# {time.asctime()} - Reading files.")
    parts = list()
    for part in range(1, n_parts + 1):
        print(f"# {time.asctime()} - Reading part {part}.")
        input_file = os.path.join(f"data/{sample_name}/{sample_name}.{part}.STAR.filtered2.csv.gz")

        d = pd.read_csv(
            input_file,
            header=None, sep=',', error_bad_lines=False, engine='c',
            names=['read', 'r1', 'r2', 'umi', 'gene', 'pos'], compression='gzip')
        parts.append(d)
        print(f"# {time.asctime()} - Done with part {part}. {d.shape[0]} lines.")

    return pd.concat(parts)


def to_pickle(obj, name, array=True):
    print(f"# {time.asctime()} - Saving {name} to pickle.")
    if array:
        pickle.dump(obj.values, open(input_prefix + f"{name}.values.pickle", 'wb'), protocol=-1)
    pickle.dump(obj, open(input_prefix + f"{name}.pickle", 'wb'), protocol=-1)


def from_pickle(key, array=False):
    print(f"# {time.asctime()} - Loading {key} from pickle.")
    if array:
        return pickle.load(open(input_prefix + f"{key}.values.pickle", 'rb'))
    return pickle.load(open(input_prefix + f"{key}.pickle", 'rb'))

    # df = pickle.load(open(input_prefix + "all.pickle", 'rb'))


def add_colorbar_to_axis(collection, label=None, position="right", size="5%", pad=0.05):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(collection.axes)
    cax = divider.append_axes(position, size=size, pad=pad)
    cbar = plt.colorbar(mappable=collection, cax=cax, label=label, alpha=1)
    cbar.solids.set_edgecolor("face")
    cbar.solids.set_rasterized(True)


def f(x, m, b):
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
