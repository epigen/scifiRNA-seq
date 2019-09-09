#! /bin/env python

import sys
import os
import time
from argparse import ArgumentParser

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc

from ngs_toolkit.general import query_biomart


def parse_args():
    parser = ArgumentParser()
    parser.add_argument(dest="input_h5ad")
    parser.add_argument(dest="output_prefix")
    parser.add_argument("--name", dest="name")
    parser.add_argument("--r1-attributes", dest="r1_attributes", type=str)
    parser.add_argument("--droplet-column", dest="droplet_column", default="r2")
    parser.add_argument("--well-column", dest="well_column", default="well")
    parser.add_argument("--species-mixture", dest="species_mixture", action="store_true")

    # Example:
    # args = parser.parse_args("--name SCI024_Tcell_s.exon.20190617. results/SCI024_Tcell_s.exon.h5ad results/SCI024_Tcell_s.exon.20190617. --r1-attributes donor_id".split(" "))
    # args = parser.parse_args("--name PD190_sixlines. results/PD190_sixlines.h5ad results/PD190_sixlines. --species-mixture --r1-attributes plate_well,cell_line".split(" "))
    args = parser.parse_args()

    if args.name is None:
        args.name = args.output_prefix
    args.r1_attributes = args.r1_attributes.split(",")
    return args


def main():
    global args
    args = parse_args()

    # barcode annotations
    annotation_file = os.path.join("metadata", "sciRNA-seq.PD190_sixlines.oligos_2019-09-05.csv")
    annotation = pd.read_csv(annotation_file)

    # convenience
    gene_set_libraries = [
        'Human_Gene_Atlas', 'ARCHS4_Tissues', 'WikiPathways_2019_Human',
        'NCI-Nature_2016', 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X', 'GO_Biological_Process_2018']

    # read h5ad file
    sc.settings.n_jobs = -1
    sc.settings.figdir = os.path.dirname(args.output_prefix)
    sc.settings.set_figure_params(dpi=300, dpi_save=300, format='svg')
    print(f"# {time.asctime()} - Reading input data.")
    adata = sc.read(args.input_h5ad, cache=True)

    # Annotate with gene names instead of Ensembl IDs
    print(f"# {time.asctime()} - Annotating genes.")
    if args.species_mixture:
        adata.var.loc[:, "species"] = pd.Series(
            adata.var_names.str.startswith("ENSMUS"),
            index=adata.var.index
        ).replace(True, "mouse").replace(False, "human")
    human_m = query_biomart(attributes=["ensembl_gene_id", "external_gene_name"], species="hsapiens", ensembl_version='grch38')
    v = adata.var.join(human_m.set_index("ensembl_gene_id"))

    if args.species_mixture:
        mouse_m = query_biomart(attributes=["ensembl_gene_id", "external_gene_name"], species="mmusculus", ensembl_version='grcm38')
        v.update(mouse_m.set_index("ensembl_gene_id"))
    adata.var.index = v['external_gene_name'].fillna(v.index.to_series()).values
    adata.var_names_make_unique()

    # QC
    # sc.pl.highest_expr_genes(a, n_top=20)
    adata.var.loc[:, 'mito'] = adata.var_names.str.contains(r'^MT-', case=False)
    adata.obs.loc[:, 'percent_mito'] = np.sum(
        adata[:, adata.var['mito']].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    adata.var.loc[:, 'ribo'] = adata.var_names.str.contains(r'^RP', case=False)
    adata.obs.loc[:, 'percent_ribo'] = np.sum(
        adata[:, adata.var['ribo']].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    adata.obs.loc[:, 'n_counts'] = adata.X.sum(axis=1).A1
    adata.obs.loc[:, 'log_counts'] = np.log10(adata.obs.loc[:, 'n_counts'])
    adata.obs.loc[:, 'n_genes'] = (adata.X != 0).sum(1).A1

    # Filter
    print(f"# {time.asctime()} - Filtering.")
    sc.pp.filter_cells(
        adata, min_counts=100)
    sc.pp.filter_cells(
        adata, max_counts=8000)
    sc.pp.filter_genes(adata, min_counts=20)
    # sc.pp.filter_genes(adata, max_counts=5000)
    print(f"Kept {adata.shape[0]} cells and {adata.shape[1]} genes.")

    # # remove cells with extreme mitochondial/ribosomal expression

    # # visualize
    # sc.pl.violin(adata, ['n_genes', 'n_counts', 'percent_mito', 'percent_ribo'], jitter=0.4, multi_panel=True)

    # Add experiment-specific variables
    print(f"# {time.asctime()} - Adding experiment-specific variables.")
    # info = adata.obs.index.to_series().str.split("-").apply(pd.Series)
    # info.columns = ['plate', 'well', 'droplet']
    info = adata.obs.index.str.slice(0, 3)
    adata.obs = adata.obs.assign(plate_well=info)
    # remove cells not matching annotation
    if adata.obs['plate_well'].isnull().sum() > 0:
        print(f"# {time.asctime()} - Warning: not all cells matched plate_well annotation.")
        adata = adata[~adata.obs['plate_well'].isnull(), :]
    adata.obs = adata.obs.merge(annotation[args.r1_attributes], on=["plate_well"], validate='many_to_one').set_index(adata.obs.index)

    if args.r1_attributes == ['plate_well', 'cell_line']:
        adata.obs = adata.obs.assign(
            species=(adata.obs['cell_line'] == "3T3").replace(True, "mouse").replace(False, "human"))

    adata.X = adata.X.astype(np.float)
    adata.raw = adata
    sc.write(args.input_h5ad.replace(".h5ad", ".filtered.h5ad"), adata)
    adata = sc.read(args.input_h5ad.replace(".h5ad", ".filtered.h5ad"), cache=True)

    a = adata.copy()
    # gene_count = pd.Series(a.X.sum(0).A1, index=a.var.index).sort_values()

    # Normalize
    sc.pp.normalize_per_cell(a)
    sc.pp.log1p(a)

    sc.tl.rank_genes_groups(a, 'cell_line', method='t-test_overestim_var', n_genes=50, use_raw=False)
    result = a.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    diff = pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals', 'logfoldchanges', 'scores']})
    diff.to_csv(args.output_prefix + "cell_line.cluster_comparison.top_values.csv", index=False)

    from ngs_toolkit.general import enrichr

    res = list()
    for cell_line in diff.columns[diff.columns.str.endswith("_n")]:
        res.append(
            enrichr(
                diff.rename(columns={cell_line: "gene_name"}),
                gene_set_libraries=['ARCHS4_Cell-lines']).assign(cell_line=cell_line))
    res = pd.concat(res)
    g = res.set_index("description").groupby(['cell_line'])['combined_score'].nlargest(5)
    print(g)

    # Reduce variables
    sc.pp.highly_variable_genes(a, flavor="seurat", min_disp=1, max_mean=1)  # , n_top_genes=100
    # sc.pp.highly_variable_genes(a, flavor="cell_ranger", n_top_genes=100, batch_key="species")  # , n_top_genes=100
    sc.pl.highly_variable_genes(a)

    # sc.pp.scale(a)
    # sc.pp.highly_variable_genes(a, flavor="seurat", min_disp=1.5, max_disp=1e5, min_mean=1, max_mean=1e5, n_bins=20)
    # sc.pl.highly_variable_genes(a)

    print(f"Found {a.var.highly_variable.sum()} highly variable genes.")

    sc.pp.scale(a)
    sc.pp.pca(a, svd_solver='arpack', zero_center=None, use_highly_variable=True)

    sc.pl.pca_variance_ratio(a, log=True)

    # Manifold
    # sc.pp.neighbors(a, use_rep="X_pca", n_neighbors=20, metric="correlation")
    sc.pp.neighbors(a, use_rep="X_pca", n_neighbors=20)
    sc.tl.umap(a)
    sc.tl.diffmap(a)

    # Cluster
    sc.tl.leiden(a, resolution=0.3)
    sc.pl.umap(
        a, color=['leiden', 'log_counts'], palette='tab20c', save=args.name + "leiden.cells_per_cluster.svg")
    sc.pl.diffmap(
        a, color=['leiden', 'log_counts'], palette='tab20c', save=args.name + "leiden.cells_per_cluster.svg")

    a.uns['iroot'] = np.argmin(a.obsm['X_diffmap'][0])
    sc.tl.dpt(a, n_branchings=2)
    sc.pl.dpt_groups_pseudotime(a)
    sc.pl.diffmap(a, color=['leiden', 'dpt_pseudotime'])

    # Differential genes
    sc.tl.rank_genes_groups(a, 'leiden', method='t-test', n_genes=1e6, use_raw=False)
    sc.pl.rank_genes_groups(a, n_genes=25, sharey=False, save=args.name + "leiden.svg")

    result = a.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    diff = pd.DataFrame(
        {group + '_' + key[:1]: result[key][group]
            for group in groups for key in ['names', 'pvals', 'logfoldchanges', 'scores']})
    diff.to_csv(args.output_prefix + "leiden.cluster_comparison.top_values.csv", index=False)

    # Enrichment analysis of clusters
    import gseapy as gp
    from ngs_toolkit.analysis import Analysis

    enrichr = list()
    for i in sorted(set(a.obs['leiden'].astype(int))):
        print(i)
        enrichr.append(gp.enrichr(
            gene_list=diff[f"{i}_n"].head(200),
            gene_sets=gene_set_libraries,
            cutoff=0.5).results.assign(comparison_name=i))
    enrichr = pd.concat(enrichr)
    enrichr = enrichr.rename(columns={
        "P-value": "p_value", "Term": 'description',
        "Gene_set": 'gene_set_library', 'Combined Score': 'combined_score',
        'Z-score': 'z_score'})
    enrichr.to_csv(args.output_prefix + "leiden.cluster_enrichments.csv", index=False)
    n = Analysis(genome='hg38')
    n.enrichment_results = {"enrichr": enrichr}
    n.plot_differential_enrichment(
        steps=['enrichr'],
        plot_types=['heatmap'],
        output_dir=".",
        output_prefix=args.output_prefix + "leiden", top_n=30)
    a.uns['enrichr'] = enrichr

    # Plot some top genes
    genes = diff.loc[:, diff.columns.str.endswith("_n")].head().T.stack().values
    sc.pl.dotplot(a, var_names=genes, groupby='leiden', use_raw=False, save=args.name + ".leiden.svg")
    sc.pl.stacked_violin(a, var_names=genes, groupby='leiden', use_raw=False, log=True, save=args.name + ".leiden.svg")
    genes = diff.loc[:, diff.columns.str.endswith("_n")].head(20).T.stack().values
    sc.pl.matrixplot(a, var_names=genes, groupby='leiden', use_raw=False, save=args.name + ".leiden.svg")

    # Save processed data
    sc.write(args.name + "processed.h5ad", a)
    a = sc.read(args.name + "processed.h5ad", cache=True)

    # # Donor-specific analysis
    sc.tl.rank_genes_groups(a, groupby='sex', method='t-test', n_genes=1e6, use_raw=False)
    sc.pl.rank_genes_groups(a, n_genes=25, sharey=False, save=args.name + "sex.svg")

    sc.tl.rank_genes_groups(a, groupby='donor_id', method='t-test', n_genes=1e6, use_raw=False)
    sc.pl.rank_genes_groups(a, n_genes=25, sharey=False, save=args.name + "donor_id.svg")

    # sc.pl.umap(a, color=['leiden', 'sex', 'donor_id'])

    # # sex remaining donors
    # url = "https://raw.githubusercontent.com/broadinstitute/inferCNV_examples/master/__gene_position_data/gencode_v19_gene_pos.txt"
    # gene_order = pd.read_csv(url, sep="\t", header=None)
    # gene_order.columns = ['gene', 'chr', 'start', 'end']

    # for chrom in ['X', 'Y']:
    #     gene_s = gene_order.loc[gene_order['chr'].str.contains(chrom)]
    #     g = gene_s.loc[gene_s['gene'].isin(adata.var.index.tolist()), 'gene'].squeeze().tolist()
    #     q = adata[:, g].X
    #     if isinstance(q, scipy.sparse.csr_matrix):
    #         q = q.todense()
    #     adata.obs.loc[:, chrom + "_expression"] = pd.DataFrame(q, index=adata.obs.index, columns=g).T.sum(0)
    # adata.obs.loc[:, 'sex_ratio'] = np.log2(adata.obs['Y_expression'] / adata.obs['X_expression'])
    # adata.obs.loc[:, 'predicted_sex'] = (adata.obs['sex_ratio'] >= 0).replace(True, "Male").replace(False, "Female")

    # from sklearn.ensemble import RandomForestClassifier
    # from sklearn.model_selection import cross_val_score
    # from sklearn.model_selection import KFold
    # model = RandomForestClassifier(n_estimators=100)
    # train = adata[~adata.obs.sex.isnull(), :]
    # X = train.X
    # Y = train.obs.sex.replace("Male", 0).replace("Female", 1)
    # model.fit(X, Y)
    # a.obs.loc[:, 'predicted_sex'] = pd.Series(model.predict(a.X), index=a.obs.index).replace(0, "Male").replace(1, "Female")
    # scores = cross_val_score(model, X, Y, cv=5)

    # preds = list()
    # kf = KFold(n_splits=5)
    # for i, (train, test) in enumerate(kf.split(X, Y)):
    #     print(i)
    #     model.fit(X[train, :], Y[train])
    #     preds.append(pd.Series(model.predict(X[test, :]), index=adata.obs.index[test]))
    # preds = pd.DataFrame(preds).sum(0)
    # a.obs.loc[:, 'predicted_sex_cv'] = preds.replace(0, "Male").replace(1, "Female")

    # sc.pl.umap(a, color=['log_counts', 'donor_id', 'sex', 'predicted_sex'])

    # # Try to impute
    # # Denoise with ZINB
    # a2 = adata.copy()
    # from dca.api import dca
    # model = dca(
    #     a2,
    #     mode="denoise", ae_type="nb",
    #     return_info=True, return_model=True)
    # sc.pp.normalize_per_cell(a2)
    # sc.pp.log1p(a2)
    # sc.pp.highly_variable_genes(a2, flavor="cell_ranger", n_top_genes=500)
    # a2.raw = a2
    # sc.pp.scale(a2)
    # sc.pp.pca(a2, use_highly_variable=True)
    # sc.pp.neighbors(a2)
    # sc.tl.umap(a2)

    # # Plot markers
    # sc.pl.umap(a, color=['donor_id', 'log_counts'], save=output_prefix + "umap.svg")

    # markers = set(["CD3E", "CD3D", "CD3G", "CD4", "CD8A", "NKG7", "IL32", "GZMB", "NCR1", "IFNG", "ITGA4", "LCK"] +
    #               ['CD8A', 'CD4', 'NKG7', "ZAP70", "LCK", ] +
    #               ["CD36", "CD27", "CD37", "CD38"] +
    #               ["FOXP1", "FOS", "JUN", "JUNB", "BCL6", "RUNX1", "IRF4" "BATF", "STAT3", "ATF2", "CEBPB", "FASL", "PRDM1"]
    #               +
    #               ['IL7R', 'CD79A', 'MS4A1', 'CD8A', 'CD8B', 'LYZ', 'CD14', 'EPO',
    #                'LGALS3', 'S100A8', 'GNLY', 'NKG7', 'KLRB1',
    #                'FCGR3A', 'MS4A7', 'FCER1A', 'CST3', 'PPBP',
    #                'PTPRC']
    #               )
    # available = [x for x in markers if x in a.var_names]
    # # sc.pl.pca(a, color=['n_counts'] + available)
    # sc.pl.umap(a, color=['donor_id', 'log_counts'] + available, save=output_prefix + "umap.markers.svg")

    # # sc.pl.pca(a, color=['leiden', 'donor_id', 'log_counts'])
    # sc.pl.umap(a, color=['leiden', 'donor_id', 'log_counts'], save=output_prefix + "umap.leiden.svg")
    # sc.pl.umap(a, color=['leiden', 'donor_id', 'log_counts'] + available, save=output_prefix + "umap.leiden_markers.svg")


if __name__ == "__main__":
    sns.set(context="paper", style="ticks", palette="colorblind", color_codes=True)
    matplotlib.rcParams["svg.fonttype"] = "none"
    # Don't use LaTeX for rendering
    matplotlib.rcParams["text.usetex"] = False
    sys.exit(main())
