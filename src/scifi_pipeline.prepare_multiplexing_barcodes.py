#!/usr/bin/env python

import os
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import yaml

from sklearn import preprocessing


sns.set(context="paper", style="ticks", palette="colorblind", color_codes=True)
matplotlib.rcParams["svg.fonttype"] = "none"
# Don't use LaTeX for rendering
matplotlib.rcParams["text.usetex"] = False
matplotlib.use("Agg")


SEQ_COLS = ["fixed_base1_N", "index_seq", "fixed_base2_V"]


def main():
    config = yaml.safe_load(
        open(os.path.join("metadata", "multiplexing_specification.yaml"), 'r'))

    for run_name, run_config in config.items():

        annotations = list()
        for sample_name, sample in run_config.items():
            annotations.append(annotate(sample_name, sample))
        annot = pd.concat(annotations, axis=0)

        annot[["sample_name", "multiplexing_barcode", "combinatorial_barcode"]].to_csv(
            os.path.join("metadata", run_name + "_multiplexing_barcodes.tsv"),
            sep="\t", index=False
        )

        annot.set_index("sample_name").to_csv(
            os.path.join(
                "metadata",
                f"sciRNA-seq.{run_name}.oligos_2019-10-31.full_annotation.csv")
        )

        for (sample_name, sample), annot2 in zip(run_config.items(), annotations):
            # Remove the multiplexing barcode from the name
            attrs = pd.Series(
                ['combinatorial_barcode'] + [sample['well_col'][0]] + list(sample['attributes'].keys()))
            pretty_attrs = pd.Series(
                ['combinatorial_barcode'] + [sample['well_col'][1]] + list(sample['attributes'].keys())).str.lower()

            c = annot2.drop(["multiplexing_number", "multiplexing_barcode"], axis=1)
            c['sample_name'] = c['sample_name'].str.replace(r"_\d\d$", "")
            c = c.drop_duplicates().set_index("sample_name")[attrs]
            c.columns = pretty_attrs
            date = "2019-10-31"
            prefix = os.path.join("metadata", f"sciRNA-seq.{sample_name}.oligos_{date}")
            c.to_csv(prefix + ".csv")
            if sample['attributes'].values():
                plot_plate(c, attrs=sample['attributes'].values(), output_prefix=prefix)


def annotate(sample_name, sample):
    d = pd.read_excel(
        os.path.join(
            "metadata", "original", sample['annotation']
        )
    )
    for col in d.columns:
        d.loc[:, col] = d[col].astype(str).str.strip().str.replace("-", "")
    d.loc[:, "combinatorial_barcode"] = d[SEQ_COLS].sum(axis=1)
    seqs = [(b, str(i).zfill(2)) for i, b in enumerate(sample['i7_seqs'], start=1)]
    d = pd.concat([
        d.assign(multiplexing_barcode=b, multiplexing_number=i)
        for b, i in seqs]
    ).reset_index(drop=True)
    annot = d.assign(
        sample_name=sample_name
        + "_"
        + d[sample['well_col'][0]].str.strip()
        + "_"
        + d["multiplexing_number"]
    )
    for k, v in sample['attributes'].items():
        annot[v] = d[k]
    return annot


def plot_plate(annot, attrs, output_prefix, annotate=False):
    def join(x):
        # from collections import OrderedDict
        return ' '.join(np.unique(x))

    annot['row'] = annot['plate_well'].str.slice(0, 1)
    annot['col'] = annot['plate_well'].str.slice(1)
    annot['label'] = annot[attrs].apply(join, axis=1)
    counts = annot['label'].value_counts()

    le = preprocessing.LabelEncoder()
    le.fit(annot['label'])

    plate = annot.pivot_table(
        index='row', columns='col', values='label', aggfunc=join)

    plate_num = pd.DataFrame(
        le.transform(plate.values.flatten()).reshape(plate.shape),
        index=plate.index, columns=plate.columns)

    kwargs = dict()
    if annotate:
        kwargs.update({"annot": plate.values, "fmt": ""})
    fig, axis = plt.subplots(1, 1, figsize=(12, 8))
    sns.heatmap(
        plate_num, cmap="tab20c", square=True, cbar=False, ax=axis, **kwargs)
    fig.savefig(output_prefix + ".plate_design.svg", dpi=300, bbox_inches="tight")

    fig, axis = plt.subplots(1, 1, figsize=(6, 3))
    sns.barplot(counts, counts.index, orient="horiz", ax=axis)
    fig.savefig(output_prefix + ".group_counts.svg", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        sys.exit(1)
