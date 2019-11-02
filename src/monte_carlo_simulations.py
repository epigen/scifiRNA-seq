#! /bin/env python

import os
import itertools

import numpy as np
import pandas as pd
import scipy

from joblib import Parallel, delayed
from tqdm import tqdm

import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns


sns.set(context="paper", style="ticks", palette="colorblind", color_codes=True)
matplotlib.rcParams["svg.fonttype"] = "none"
# Don't use LaTeX for rendering
matplotlib.rcParams["text.usetex"] = False


def get_collision_rate(r1, r2, n, psi=1.0):
    x = np.zeros(shape=(n, 2), dtype=float)
    # sci-RNA barcode
    x[:, 0] = np.random.randint(low=0, high=r1, size=n)
    # droplet barcode
    x[:, 1] = np.random.randint(low=0, high=r2, size=n)
    # simulate droplet loading
    p = scipy.stats.bernoulli(psi).rvs(n)
    x[p != 1, 1] = np.nan

    # count collisions and get rate
    xx = pd.DataFrame(x, columns=['r1', 'r2'])
    xc = xx.groupby('r2')['r1'].value_counts()
    return (xc != 1).sum() / xc.shape[0]


def get_collision_rate_numpy(r1, r2, n, psi=1.0):
    x = np.zeros(shape=(n, 2), dtype=float)
    # sci-RNA barcode
    x[:, 0] = np.random.randint(low=0, high=r1, size=n)
    # droplet barcode
    x[:, 1] = np.random.randint(low=0, high=r2, size=n)
    # # simulate droplet loading
    p = scipy.stats.bernoulli(psi).rvs(n)
    x[p != 1, 0] = np.nan
    # # simulate empty droplet
    p = scipy.stats.bernoulli(psi).rvs(n)
    x[p != 1, 1] = np.nan

    # count collisions and get rate
    u, c = np.unique(x, axis=0, return_counts=True)
    return (c > 1).sum() / n


def approximate_collision_rate(r1, r2, n):
    return (r2 / n) * r1


# Monte Carlo simulations of scifi-RNA-seq method
iterations = 100

loading_concentrations = np.array([1000, 10000,  125000,  250000,  500e3, 1e6, 3e6], dtype=float)
loading_concentrations *= 1.53
loading_concentrations = loading_concentrations.astype(int)
r1_barcodes = [1, 96, 96 * 4, 96 * 16]
r2_barcodes = [737000, int(3e6)]
psi_range = [0.1]

# # simulate combinatorial indexing

results = Parallel(n_jobs=-1)(
    delayed(get_collision_rate)(r1, r2, n)
    for r1 in tqdm(r1_barcodes)
    for r2 in r2_barcodes
    for n in loading_concentrations
    for psi in psi_range
    for _ in range(iterations))

df = pd.DataFrame(
    list(itertools.product(r1_barcodes, r2_barcodes, loading_concentrations, psi_range, range(iterations))),
    columns=['r1_barcodes', 'r2_barcodes', 'loaded_nuclei', 'psi', 'iteration'])
df.loc[:, 'collision_rate'] = results

# df.to_csv(os.path.join("results", "monte_carlo_simulations", "monte_carlo_simulations.csv"), index=False)
df = pd.read_csv(os.path.join("results", "monte_carlo_simulations", "monte_carlo_simulations.csv"))


df_agg = (
    df.groupby(['r1_barcodes', 'r2_barcodes', 'loaded_nuclei', 'psi'])
    ['collision_rate'].agg([max, np.mean, np.median, np.std, np.var]).reset_index())  # .apply(pd.Series.describe))


cols = len(r2_barcodes)
fig, axis = plt.subplots(
    2, cols, figsize=(cols * 4, 2 * 4),
    squeeze=False, sharey="row", sharex="row")
for i, r2 in enumerate(r2_barcodes):
    for r1 in r1_barcodes:
        d = df_agg.loc[(df_agg['r1_barcodes'] == r1) & (df_agg['r2_barcodes'] == r2)]
        for ax in axis[:, i]:
            ax.plot(
                d['loaded_nuclei'],
                d['mean'],
                linestyle="--",
                marker="o",
                label=f"{r1} round1 barcodes")
            ax.fill_between(
                d['loaded_nuclei'],
                d['mean'] - d['std'],
                d['mean'] + d['std'],
                alpha=0.2)
for ax in axis.flatten():
    ax.set_xlabel("Nuclei loaded")
    ax.set_ylabel("Collision probability")
    ax.axhline(0.1, linestyle="--", color="grey")
for i, ax in enumerate(axis[0, :]):
    ax.legend()
    ax.set_title(f"{r2_barcodes[i]} round2 barcodes")
for i, ax in enumerate(axis[1, :]):
    ax.loglog()
    ax.set_title(f"{r2_barcodes[i]} round2 barcodes")

fig.savefig(
    os.path.join("results", "monte_carlo_simulations", "monte_carlo_simulations.svg"),
    bbox_inches="tight", tight_layout=True)
fig.savefig(
    os.path.join("results", "monte_carlo_simulations", "monte_carlo_simulations.png"),
    bbox_inches="tight", tight_layout=True, dpi=300)


# Compare prediction with observed from 1 round1 barcode (10X data)
droplet_counts = pd.read_csv(os.path.join("metadata", "droplet_counts.csv"))
droplet_counts = droplet_counts.rename(columns={"count": "Droplets"})
observed_collisions = droplet_counts.groupby("loaded_nuclei").apply(
    lambda x: x.loc[(x['cells_per_droplet'] > 1), 'Droplets'].sum() / x.loc[:, 'Droplets'].sum())

p = df_agg.loc[(df_agg['r1_barcodes'] == 1) & (df_agg['r2_barcodes'] == r2_barcodes[0])].set_index("loaded_nuclei")['mean']
p.name = "predicted"
joint = observed_collisions.to_frame(name="observed").join(p).dropna()

fig, axis = plt.subplots(
    1, 2, figsize=(2 * 4, 1 * 4))
for ax in axis:
    ax.plot(
        joint['predicted'],
        joint['observed'],
        linestyle="--",
        marker="o")
    ax.plot(
        (0, 1),
        (0, 1),
        linestyle="--",
        color="grey")
    ax.set_xlabel("Predicted")
    ax.set_ylabel("Observed")
axis[1].loglog()

fig.savefig(
    os.path.join("results", "monte_carlo_simulations", "monte_carlo_simulations.comparison_to_observed.svg"),
    bbox_inches="tight", tight_layout=True)
fig.savefig(
    os.path.join("results", "monte_carlo_simulations", "monte_carlo_simulations.comparison_to_observed.png"),
    bbox_inches="tight", tight_layout=True, dpi=300)
