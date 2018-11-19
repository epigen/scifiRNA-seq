#!/usr/bin/env python

"""
"""


import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import string


pd.set_option("date_dayfirst", True)
sns.set(context="paper", style="white", palette="pastel", color_codes=True)
sns.set_palette(sns.color_palette("colorblind"))
matplotlib.rcParams["svg.fonttype"] = "none"
matplotlib.rc('text', usetex=False)


# params
df = pd.read_csv(os.path.join("metadata", "annotation.csv"))
df = df[df['experimental_batch'].isin(["SCI012", "SCI013"])]
cell_barcodes = ['round1', 'round2']
mismatches = 0

annot_12 = pd.read_csv(os.path.join("metadata", "/scratch/lab_bock/shared/projects/sci-rna/metadata/sciRNA-seq.SCI012.oligos_2018-11-14.csv"))
annot_12['plate_half'] = annot_12['plate_well'].str.contains("E|F|G|H").astype(int)
annot_12.loc[annot_12['plate_well'].str.contains("A|B"), 'plate_quarter'] = 0
annot_12.loc[annot_12['plate_well'].str.contains("C|D"), 'plate_quarter'] = 1
annot_12.loc[annot_12['plate_well'].str.contains("E|F"), 'plate_quarter'] = 2
annot_12.loc[annot_12['plate_well'].str.contains("G|H"), 'plate_quarter'] = 3
annot_12['plate_quarter'] = annot_12['plate_quarter'].astype(int)
annot_12['plate_octo'] = np.nan
q = -1
for l in string.ascii_uppercase[:8]:
    for n in [str(x) for x in range(1, 13)]:
        if n in ['1', '7']:
            q += 1
        if len(n) == 1:
            n = '0' + n
        annot_12.loc[annot_12['plate_well'] == (l + n), 'plate_octo'] = q
annot_12['plate_octo'] = annot_12['plate_octo'].astype(int)
annot_12_1 = annot_12.loc[annot_12['barcode_type'] == "round1", ['barcode_sequence', 'plate_well', 'plate_half', 'plate_quarter', 'plate_octo']]
annot_12_1 = annot_12_1.rename(columns={"barcode_sequence": 'round1', "plate_well": "round1_well"})

annot_13 = pd.read_csv(os.path.join("metadata", "/scratch/lab_bock/shared/projects/sci-rna/metadata/sciRNA-seq.SCI013.oligos_2018-11-14.csv"))
annot_13['plate_half'] = annot_13['plate_well'].str.contains("E|F|G|H").astype(int)
annot_13.loc[annot_13['plate_well'].str.contains("A|B"), 'plate_quarter'] = 0
annot_13.loc[annot_13['plate_well'].str.contains("C|D"), 'plate_quarter'] = 1
annot_13.loc[annot_13['plate_well'].str.contains("E|F"), 'plate_quarter'] = 2
annot_13.loc[annot_13['plate_well'].str.contains("G|H"), 'plate_quarter'] = 3
annot_13['plate_quarter'] = annot_13['plate_quarter'].astype(int)
annot_13['plate_octo'] = np.nan
q = -1
for l in string.ascii_uppercase[:8]:
    for n in [str(x) for x in range(1, 13)]:
        if n in ['1', '7']:
            q += 1
        if len(n) == 1:
            n = '0' + n
        annot_13.loc[annot_13['plate_well'] == (l + n), 'plate_octo'] = q
annot_13['plate_octo'] = annot_13['plate_octo'].astype(int)
annot_13_1 = annot_13.loc[annot_13['barcode_type'] == "round1", ['barcode_sequence', 'plate_well', 'plate_half', 'plate_quarter', 'plate_octo']]
annot_13_1 = annot_13_1.rename(columns={"barcode_sequence": 'round1', "plate_well": "round1_well"})


labels = {
    "sci-RNA-seq_SCI012_12-nuclei_various-cycles": {
        {k: v for k, v in zip(range(16), ["standard"] * 8 + ['cycling'] * 8)}},
    "sci-RNA-seq_SCI012_24-nuclei_various-cycles": {
        {k: v for k, v in zip(range(16), ["standard"] * 8 + ['cycling'] * 8)}},
    "sci-RNA-seq_SCI012_12-nuclei_various-enzymes": {
        {k: v for k, v in zip(range(16), ["ProtoscriptII"] * 4 + ['Maxima Hminus'] * 4 + ['Superscript IV'] * 4 + ['Revert Aid minus'] * 4)}},
    "sci-RNA-seq_SCI012_24-nuclei_various-enzymes": {
        {k: v for k, v in zip(range(16), ["ProtoscriptII"] * 4 + ['Maxima Hminus'] * 4 + ['Superscript IV'] * 4 + ['Revert Aid minus'] * 4)}},
    "sci-RNA-seq_SCI013_12-cells_various-rtprimers": {
        {k: v for k, v in zip(range(16), ["Anchored"] * 4 + ['Not anchored'] * 4 + ['NA'] * 8)}},
    "sci-RNA-seq_SCI013_12-nuclei_various-rtprimers": {
        {k: v for k, v in zip(range(16), ["Anchored"] * 4 + ['Not anchored'] * 4 + ['NA'] * 8)}},
}

# start
u_res = pd.DataFrame()
mr_res = pd.DataFrame()
dr_res = pd.DataFrame()
for sample in df['sample_name']:
    prefix = "{}.{}mis.".format(sample, mismatches)
    print(prefix)

    if "SCI012" in sample:
        a = annot_12_1
    else:
        a = annot_13_1

    # get umi counts
    umi_count_file = os.path.join("barcodes", prefix + "barcode_gene_umi_count.clean.csv")
    umi = pd.read_csv(umi_count_file)
    # # add well info
    umi = umi.merge(a, on="round1")
    for v in range(16):
        umi.loc[umi['plate_octo'] == v, "plate_octo"] = labels[sample][v]
    umi['sample_name'] = sample
    u_res = u_res.append(umi)

    # get mapping rate per well
    mapping_rate_file = os.path.join("barcodes", prefix + "mapping_rate.per_well.csv")
    mr = pd.read_csv(mapping_rate_file, index_col=0)
    mr = pd.merge(mr.reset_index(), a, on=["round1_well"])
    for v in range(16):
        mr.loc[mr['plate_octo'] == v, "plate_octo"] = labels[sample][v]
    mr['sample_name'] = sample
    mr_res = mr_res.append(mr)

    dup_rate_file = os.path.join("barcodes", prefix + "barcode_umi_dups.per_well.csv")
    dr = pd.read_csv(dup_rate_file, index_col=0)
    dr = pd.merge(dr.reset_index(), a, on=["round1_well"])
    for v in range(16):
        dr.loc[dr['plate_octo'] == v, "plate_octo"] = labels[sample][v]
    dr['sample_name'] = sample
    dr_res = dr_res.append(dr)


# # umis per cell / well
umi_cell = u_res.groupby(['sample_name'] + cell_barcodes + ['plate_octo'])['umi'].sum().reset_index()
umi_cell_filtered = umi_cell[umi_cell['umi'] >= 500]

fig, axis = plt.subplots(2, 1, figsize=(6, 4 * 2))
axis[0].set_title("All barcodes")
sns.barplot(
    data=umi_cell, y="sample_name", x="umi",
    hue="plate_octo", orient="horiz", ax=axis[0], hue_order=umi_cell['plate_octo'].unique(), capsize=.1)
axis[1].set_title("Real cells")
sns.barplot(
    data=umi_cell_filtered, y="sample_name", x="umi",
    hue="plate_octo", orient="horiz", ax=axis[1], hue_order=umi_cell['plate_octo'].unique(), capsize=.1)
for ax in axis:
    ax.set_xlabel("UMIs per cell")
    ax.set_ylabel("Experiment")
sns.despine(fig)
fig.savefig(os.path.join("results", "SCI012-13.comparisons.UMIs_per_cell.svg"), dpi=300, bbox_inches="tight", tight_layout=True)

# # genes per cell / well
gene_cell = u_res.groupby(['sample_name'] + cell_barcodes + ['plate_octo'])['gene'].nunique().reset_index()
gene_cell_filtered = gene_cell.loc[umi_cell_filtered.index]

fig, axis = plt.subplots(2, 1, figsize=(6, 4 * 2))
axis[0].set_title("All barcodes")
sns.barplot(
    data=gene_cell, y="sample_name", x="gene",
    hue="plate_octo", orient="horiz", ax=axis[0], hue_order=umi_cell['plate_octo'].unique(), capsize=.1)
axis[1].set_title("Real cells")
sns.barplot(
    data=gene_cell_filtered, y="sample_name", x="gene",
    hue="plate_octo", orient="horiz", ax=axis[1], hue_order=umi_cell['plate_octo'].unique(), capsize=.1)
for ax in axis:
    ax.set_xlabel("Genes per cell")
    ax.set_ylabel("Experiment")
sns.despine(fig)
fig.savefig(os.path.join("results", "SCI012-13.comparisons.genes_per_cell.svg"), dpi=300, bbox_inches="tight", tight_layout=True)


# # mapping rate per cell / well
fig, axis = plt.subplots(1, 1, figsize=(6, 4 * 1))
axis.set_title("All barcodes")
sns.barplot(
    data=mr_res, y="sample_name", x="mapping_rate",
    hue="plate_octo", orient="horiz", ax=axis, hue_order=umi_cell['plate_octo'].unique(), capsize=.1)
axis.set_xlabel("Mapping rate")
axis.set_ylabel("Experiment")
sns.despine(fig)
fig.savefig(os.path.join("results", "SCI012-13.comparisons.mapping_rate.svg"), dpi=300, bbox_inches="tight", tight_layout=True)


# # duplication rate per cell / well
fig, axis = plt.subplots(1, 1, figsize=(6, 4 * 1))
axis.set_title("All barcodes")
sns.barplot(
    data=dr_res, y="sample_name", x="unique_rate",
    hue="plate_octo", orient="horiz", ax=axis, hue_order=umi_cell['plate_octo'].unique(), capsize=.1)
axis.set_xlabel("Unique UMI rate")
axis.set_ylabel("Experiment")
sns.despine(fig)
fig.savefig(os.path.join("results", "SCI012-13.comparisons.unique_rate.svg"), dpi=300, bbox_inches="tight", tight_layout=True)


# # grnas captured per cell / well
u_res['grna'] = u_res['gene'].str.contains("Tcrlibrary|CTRL0").astype(int)
grna_cell = u_res.groupby(['sample_name'] + cell_barcodes + ['plate_octo'])['grna'].sum().reset_index()
grna_cell_filtered = u_res.loc[umi_cell_filtered.index].groupby(['sample_name'] + cell_barcodes + ['plate_octo'])['grna'].sum().reset_index()

fig, axis = plt.subplots(2, 1, figsize=(6, 4 * 2))
axis[0].set_title("All barcodes")
sns.barplot(
    data=grna_cell, y="sample_name", x="grna",
    hue="plate_octo", orient="horiz", ax=axis[0], hue_order=grna_cell['plate_octo'].unique(), capsize=.1)
axis[1].set_title("Real cells")
sns.barplot(
    data=grna_cell_filtered, y="sample_name", x="grna",
    hue="plate_octo", orient="horiz", ax=axis[1], hue_order=grna_cell['plate_octo'].unique(), capsize=.1)
for ax in axis:
    ax.set_xlabel("gRNAs UMIs per cell")
    ax.set_ylabel("Experiment")
sns.despine(fig)
fig.savefig(os.path.join("results", "SCI012-13.comparisons.gRNAs_per_cell.svg"), dpi=300, bbox_inches="tight", tight_layout=True)


# # Cas9 per cell / well
u_res['cas9'] = u_res['gene'].str.contains("Cas9").astype(int)
cas9_cell = u_res.groupby(['sample_name'] + cell_barcodes + ['plate_octo'])['cas9'].sum().reset_index()
cas9_cell_filtered = u_res.loc[umi_cell_filtered.index].groupby(['sample_name'] + cell_barcodes + ['plate_octo'])['cas9'].sum().reset_index()

fig, axis = plt.subplots(2, 1, figsize=(6, 4 * 2))
axis[0].set_title("All barcodes")
sns.barplot(
    data=cas9_cell, y="sample_name", x="cas9",
    hue="plate_octo", orient="horiz", ax=axis[0], hue_order=cas9_cell['plate_octo'].unique(), capsize=.1)
axis[1].set_title("Real cells")
sns.barplot(
    data=cas9_cell_filtered, y="sample_name", x="cas9",
    hue="plate_octo", orient="horiz", ax=axis[1], hue_order=cas9_cell['plate_octo'].unique(), capsize=.1)
for ax in axis:
    ax.set_xlabel("Cas9 UMIs per cell")
    ax.set_ylabel("Experiment")
sns.despine(fig)
fig.savefig(os.path.join("results", "SCI012-13.comparisons.Cas9_per_cell.svg"), dpi=300, bbox_inches="tight", tight_layout=True)

 # # Barcode representation

# ## molecule level
# for barcode in ['round1', 'round2']:
#     fig, axis = plt.subplots(2, len(runs), figsize=(len(runs) * 2, 2 * 2))
#     for i, mismatches in enumerate(range(2)):
#         for j, sample in enumerate(runs):
#             r = results.loc[(results['sample'] == sample) & (results['mismatches'] == mismatches), barcode].value_counts()
#             axis[i, j].plot(r.rank(), r)
#             axis[0, j].set_title(sample)
#             axis[1, j].set_xlabel("Barcode")
#         axis[i, 0].set_ylabel("Barcode abundance")
#     sns.despine(fig)
#     fig.savefig(
#         os.path.join(args.output_dir, "all_experiments.{}_barcode_abundance.lineplot.svg".format(barcode)),
#         bbox_inches="tight")

# for barcode in ['round1', 'round2']:
#     fig, axis = plt.subplots(2, 2, figsize=(2 * 4, 2 * 4))
#     for i, mismatches in enumerate(range(2)):
#         for j, sample in enumerate(runs):
#             r = results.loc[(results['sample'] == sample) & (results['mismatches'] == mismatches), barcode].value_counts()
#             axis[i, 0].plot(r.rank(), r, label=sample)
#             axis[i, 1].plot(r.rank(), r / r.sum(), label=sample)
#             for ax in axis[i, :]:
#                 ax.set_xlabel("Barcode")
#         axis[i, 0].legend()
#         axis[i, 0].set_ylabel("Barcode abundance")
#         axis[i, 1].set_ylabel("Barcode fraction")
#     sns.despine(fig)
#     fig.savefig(
#         os.path.join(args.output_dir, "all_experiments.{}_barcode_abundance.lineplot.overlay.svg".format(barcode)),
#         bbox_inches="tight")

# ## cell level
# ### rank vs abundance
# for mismatches in range(2):
#     fig, axis = plt.subplots(2, 2, figsize=(2 * 4, 2 * 4))
#     for i, barcode in enumerate(['round1', 'round2']):
#         for j, sample in enumerate(runs):
#             c = results.loc[(results['sample'] == sample) & (results['mismatches'] == mismatches), :]
#             c['cell'] = c['round1'] + c['round2']

#             # get distribution of barcodes for each cell
#             c2 = c[['round1', 'round2']].drop_duplicates().loc[:, barcode].value_counts()
#             # axis[i, 0].plot(c2.rank(), c2, label=sample)
#             axis[i, 0].plot(c2.rank(), c2 / c2.sum(), label=sample)

#             # now get same distribution for top cells only
#             top_cells = c.groupby('cell')['umi'].sum().sort_values().tail(1000).index.tolist()
#             cc = c[c['cell'].isin(top_cells)]
#             cc2 = cc[['round1', 'round2']].drop_duplicates().loc[:, barcode].value_counts()
#             # axis[i, 0].plot(cc2.rank(), cc2, label=sample)
#             axis[i, 1].plot(cc2.rank(), cc2 / cc2.sum(), label=sample)

#         axis[i, 0].legend()
#         for ax in axis[0, :]:
#             ax.set_xlabel("Round1 barcodes")
#         for ax in axis[1, :]:
#             ax.set_xlabel("Round2 barcodes")
#         for ax in axis[i, :]:
#             ax.set_ylabel("Cell fraction")
#         axis[0, 0].set_title("All cells")
#         axis[0, 1].set_title("Top 1000 cells")
#     sns.despine(fig)
#     fig.savefig(
#         os.path.join(args.output_dir, "all_experiments.mis{}.barcode_abundance.cell_level.lineplot.overlay.svg".format(mismatches)),
#         bbox_inches="tight")
# ### matched barcodes
# for mismatches in range(2):
#     fig, axis = plt.subplots(2, 2, figsize=(2 * 4, 2 * 4))
#     for i, barcode in enumerate(['round1', 'round2']):
#         bc = results[barcode].drop_duplicates().sort_values()
#         for j, sample in enumerate(runs):
#             c = results.loc[(results['sample'] == sample) & (results['mismatches'] == mismatches), :]
#             c['cell'] = c['round1'] + c['round2']

#             # get distribution of barcodes for each cell
#             c2 = c[['round1', 'round2']].drop_duplicates().loc[:, barcode].value_counts()
#             c2 = c2.loc[bc]
#             # axis[i, 0].scatter(range(c2.shape[0]), c2, label=sample)
#             axis[i, 0].scatter(range(c2.shape[0]), c2 / c2.sum(), label=sample)

#             # now get same distribution for top cells only
#             top_cells = c.groupby('cell')['umi'].sum().sort_values().tail(1000).index.tolist()
#             cc = c[c['cell'].isin(top_cells)]
#             cc2 = cc[['round1', 'round2']].drop_duplicates().loc[:, barcode].value_counts()
#             cc2 = cc2.loc[bc]
#             # axis[i, 0].scatter(range(cc2.shape[0]), cc2, label=sample)
#             axis[i, 1].scatter(range(cc2.shape[0]), cc2 / cc2.sum(), label=sample)

#         axis[i, 0].legend()
#         for ax in axis[0, :]:
#             ax.set_xlabel("Round1 barcodes")
#         for ax in axis[1, :]:
#             ax.set_xlabel("Round2 barcodes")
#         for ax in axis[i, :]:
#             ax.set_ylabel("Cell fraction")
#         axis[0, 0].set_title("All cells")
#         axis[0, 1].set_title("Top 1000 cells")
#     sns.despine(fig)
#     fig.savefig(
#         os.path.join(args.output_dir, "all_experiments.mis{}.barcode_abundance.cell_level.ordered.lineplot.overlay.svg".format(mismatches)),
#         bbox_inches="tight")

# ### Check row bias into species mixture
# annot1 = annotation.loc[annotation['barcode_type'] == 'round1', ['barcode_sequence', 'plate_well']].rename(
#         columns={"barcode_sequence": "round1", "plate_well": "round1_well"})
# annot2 = annotation.loc[annotation['barcode_type'] == 'round2', ['barcode_sequence', 'plate_well']].rename(
#         columns={"barcode_sequence": "round2", "plate_well": "round2_well"})
# for mismatches in range(2):
#     fig, axis = plt.subplots(2, len(runs), figsize=(len(runs) * 4, 2 * 4))
#     fig.suptitle("Mismatches: {}".format(mismatches))
#     for j, sample in enumerate(runs):
#         print(mismatches, sample)
#         c = results.loc[(results['sample'] == sample) & (results['mismatches'] == mismatches), :]
#         c['species'] = np.nan
#         c.loc[c['gene'].str.startswith("ENSG"), "species"] = 'human'
#         c.loc[c['gene'].str.startswith("ENSMUS"), "species"] = 'mouse'

#         sp = c.groupby(['round1', 'round2'])['species'].apply(lambda x: (x == 'human').sum() / float(x.shape[0])).reset_index()
#         sp = pd.merge(pd.merge(sp, annot1), annot2).sort_values(['round1', 'round2'])
#         sp['cell'] = sp['round1'] + sp['round2']
#         sp['round1_row'] = sp['round1_well'].str.slice(0, 1)
#         sp['round2_row'] = sp['round2_well'].str.slice(0, 1)

#         sns.violinplot(data=sp, x="species", y="round1_row", order=sp['round1_row'].drop_duplicates().sort_values(), orient="horiz", ax=axis[0, j])
#         sns.violinplot(data=sp, x="species", y="round2_row", order=sp['round2_row'].drop_duplicates().sort_values(), orient="horiz", ax=axis[1, j])

#         axis[1, j].set_xlabel("Fraction of human transcriptome")
#         axis[0, j].set_ylabel("Round1")
#         axis[1, j].set_ylabel("Round2")
#         axis[0, j].set_title(sample)

#     sns.despine(fig)
#     fig.savefig(
#         os.path.join(args.output_dir, "all_experiments.mis{}.species_mix.row_level.violinplot.svg".format(mismatches)),
#         bbox_inches="tight")


# #### top cells
# results['cell'] = results['round1'] + results['round2']
# top_cells = results.groupby(['cell'])['umi'].sum().sort_values().tail(1000).index.tolist()
# top_results = results.loc[results['cell'].isin(top_cells)]
# for mismatches in range(2):
#     fig, axis = plt.subplots(2, len(runs), figsize=(len(runs) * 4, 2 * 4))
#     fig.suptitle("Mismatches: {}".format(mismatches))
#     for j, sample in enumerate(runs):
#         print(mismatches, sample)
#         c = top_results.loc[(top_results['sample'] == sample) & (top_results['mismatches'] == mismatches), :]
#         c['species'] = np.nan
#         c.loc[c['gene'].str.startswith("ENSG"), "species"] = 'human'
#         c.loc[c['gene'].str.startswith("ENSMUS"), "species"] = 'mouse'

#         sp = c.groupby(['round1', 'round2'])['species'].apply(lambda x: (x == 'human').sum() / float(x.shape[0])).reset_index()
#         sp = pd.merge(pd.merge(sp, annot1), annot2).sort_values(['round1', 'round2'])
#         sp['cell'] = sp['round1'] + sp['round2']
#         sp['round1_row'] = sp['round1_well'].str.slice(0, 1)
#         sp['round2_row'] = sp['round2_well'].str.slice(0, 1)

#         sns.violinplot(data=sp, x="species", y="round1_row", order=sp['round1_row'].drop_duplicates().sort_values(), orient="horiz", ax=axis[0, j])
#         sns.violinplot(data=sp, x="species", y="round2_row", order=sp['round2_row'].drop_duplicates().sort_values(), orient="horiz", ax=axis[1, j])

#         axis[1, j].set_xlabel("Fraction of human transcriptome")
#         axis[0, j].set_ylabel("Round1")
#         axis[1, j].set_ylabel("Round2")
#         axis[0, j].set_title(sample)

#     sns.despine(fig)
#     fig.savefig(
#         os.path.join(args.output_dir, "all_experiments.mis{}.species_mix.row_level.violinplot.top_cells.svg".format(mismatches)),
#         bbox_inches="tight")
# # scatter
# for mismatches in range(2):
#     fig2, axis2 = plt.subplots(8, len(runs), figsize=(len(runs) * 2, 8 * 2), sharex=True, sharey=True)
#     fig2.suptitle("Mismatches: {}".format(mismatches))
#     for j, sample in enumerate(runs):
#         print(mismatches, sample)
#         c = results.loc[(results['sample'] == sample) & (results['mismatches'] == mismatches), :]
#         c['species'] = np.nan
#         c.loc[c['gene'].str.startswith("ENSG"), "species"] = 'human'
#         c.loc[c['gene'].str.startswith("ENSMUS"), "species"] = 'mouse'

#         # scatter
#         c = pd.merge(pd.merge(c, annot1), annot2).sort_values(['round1', 'round2'])
#         c['round1_row'] = c['round1_well'].str.slice(0, 1)
#         c['round2_row'] = c['round2_well'].str.slice(0, 1)
#         for i, row in enumerate(c['round1_row'].drop_duplicates().sort_values()):
#             species_count = pd.pivot_table(
#                 c[c['round1_row'] == row], index=args.cell_barcodes, columns="species", values="umi", aggfunc=sum, fill_value=0)
#             axis2[i, j].scatter(
#                 species_count.loc[:, "mouse"],
#                 species_count.loc[:, "human"], s=2, alpha=0.2, rasterized=True, label=row)
#             axis2[i, j].set_ylabel(row)
#         axis2[0, j].set_title(sample)
#     for ax in axis2.flatten():
#         ax.axhline(0, linestyle="--", color="black", alpha=0.1, zorder=100)
#         ax.axvline(0, linestyle="--", color="black", alpha=0.1, zorder=100)
#     sns.despine(fig2)
#     fig2.savefig(
#         os.path.join(args.output_dir, "all_experiments.mis{}.species_mix.row_level.scatter.svg".format(mismatches)),
#         bbox_inches="tight")
# # scatter log2
# for mismatches in range(2):
#     fig2, axis2 = plt.subplots(8, len(runs), figsize=(len(runs) * 2, 8 * 2), sharex=True, sharey=True)
#     fig2.suptitle("Mismatches: {}".format(mismatches))
#     for j, sample in enumerate(runs):
#         print(mismatches, sample)
#         c = results.loc[(results['sample'] == sample) & (results['mismatches'] == mismatches), :]
#         c['species'] = np.nan
#         c.loc[c['gene'].str.startswith("ENSG"), "species"] = 'human'
#         c.loc[c['gene'].str.startswith("ENSMUS"), "species"] = 'mouse'
#         # scatter
#         c = pd.merge(pd.merge(c, annot1), annot2).sort_values(['round1', 'round2'])
#         c['round1_row'] = c['round1_well'].str.slice(0, 1)
#         c['round2_row'] = c['round2_well'].str.slice(0, 1)
#         for i, row in enumerate(c['round1_row'].drop_duplicates().sort_values()):
#             species_count = pd.pivot_table(
#                 c[c['round1_row'] == row], index=args.cell_barcodes, columns="species", values="umi", aggfunc=sum, fill_value=0)
#             axis2[i, j].scatter(
#                 np.log2(1 + species_count.loc[:, "mouse"]),
#                 np.log2(1 + species_count.loc[:, "human"]), s=2, alpha=0.2, rasterized=True, label=row)
#             axis2[i, j].set_ylabel(row)
#         axis2[0, j].set_title(sample)
#     for ax in axis2.flatten():
#         ax.axhline(0, linestyle="--", color="black", alpha=0.1, zorder=100)
#         ax.axvline(0, linestyle="--", color="black", alpha=0.1, zorder=100)
#     sns.despine(fig2)
#     fig2.savefig(
#         os.path.join(args.output_dir, "all_experiments.mis{}.species_mix.row_level.scatter.log2.svg".format(mismatches)),
#         bbox_inches="tight")
