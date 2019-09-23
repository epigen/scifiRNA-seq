import os
import pandas as pd


# PD191_organoid_1.69M      TAAGGCGA
# PD191_organoid_1.69M      CGTACTAG
# PD192_PBMCs_1.96M         AGGCAGAA
# PD192_PBMCs_1.96M         TCCTGAGC
# PD193_5lines_beads_383k   GGACTCCT
# PD193_human/mouse_765k    CTCTCTAC
# PD193_human/mouse_1.53M   CGAGGCTG

seq_cols = ["fixed_base1_N", "index_seq", "fixed_base2_V"]

# PD191
d = pd.read_excel(
    os.path.join(
        "metadata", "original", "SCIFI-LIG-384RT-PBMCs.xlsx"
    )
)
for col in d.columns:
    d.loc[:, col] = d[col].astype(str).str.strip().str.replace("-", "")
d.loc[:, "combinatorial_barcode"] = d[seq_cols].sum(axis=1)
seqs = [("TAAGGCGA", "01"), ("CGTACTAG", "02")]
d = pd.concat(
    [d.assign(multiplexing_barcode=b, multiplexing_number=i) for b, i in seqs]
).reset_index(drop=True)
annot1 = d.assign(
    sample_name="PD191_organoid_169M"
    + "_"
    + d["384WP_well"].str.strip()
    + "_"
    + d["multiplexing_number"]
)

# PD192
d = pd.read_excel(
    os.path.join(
        "metadata", "original", "SCIFI-LIG-384RT-PBMCs.xlsx"
    )
)
for col in d.columns:
    d.loc[:, col] = d[col].astype(str).str.strip().str.replace("-", "")
d.loc[:, "combinatorial_barcode"] = d[seq_cols].sum(axis=1)
seqs = [("AGGCAGAA", "01"), ("TCCTGAGC", "02")]
d = pd.concat(
    [d.assign(multiplexing_barcode=b, multiplexing_number=i) for b, i in seqs]
).reset_index(drop=True)

annot2 = d.assign(
    sample_name="PD192_PBMCs_196M"
    + "_"
    + d["Donor_ID"]
    + "_"
    + d["384WP_well"].str.strip()
    + "_"
    + d["multiplexing_number"]
)

# PD193
d = pd.read_excel(
    os.path.join(
        "metadata", "original", "SCIFI-LIG-384RT-5lines.xlsx"
    )
)
for col in d.columns:
    d.loc[:, col] = d[col].astype(str).str.strip().str.replace("-", "")
d.loc[:, "combinatorial_barcode"] = d[seq_cols].sum(axis=1)

# Add cell line mixing experiment
annot3 = d.assign(
    sample_name="PD193_fivelines_383k"
    + "_"
    + d["Cell_line"].str.strip().str.replace("-", "")
    + "_"
    + d["384WP_well"].str.strip(),
    multiplexing_barcode="GGACTCCT"
)
annot4 = d.assign(
    sample_name="PD193_humanmouse_765k"
    + "_"
    + d["384WP_well"].str.strip(),
    multiplexing_barcode="CTCTCTAC"
)
annot5 = d.assign(
    sample_name="PD193_humanmouse_153M"
    + "_"
    + d["384WP_well"].str.strip(),
    multiplexing_barcode="CGAGGCTG"
)

annot = pd.concat([annot1, annot2, annot3, annot4, annot5], axis=0)


annot[["sample_name", "multiplexing_barcode", "combinatorial_barcode"]].to_csv(
    os.path.join("metadata", "PD191-193_multiplexing_barcodes.tsv"),
    sep="\t", index=False
)

annot.set_index("sample_name").to_csv(
    os.path.join(
        "metadata",
        "sciRNA-seq.PD191-193.oligos_2019-09-20.full_annotation.csv")
)


# Remove the multiplexing barcode from the name
c = annot1.drop(["multiplexing_number", "multiplexing_barcode"], axis=1)
c['sample_name'] = c['sample_name'].str.replace(r"_\d\d$", "")
c = c.drop_duplicates().set_index("sample_name")[['combinatorial_barcode', '384WP_well']]
c.columns = ['combinatorial_barcode', 'plate_well']
c.to_csv(
    os.path.join("metadata", "sciRNA-seq.PD191_organoid_169M.oligos_2019-09-20.csv"))

c = annot2.drop(["multiplexing_number", "multiplexing_barcode"], axis=1)
c['sample_name'] = c['sample_name'].str.replace(r"_\d\d$", "")
c = c.drop_duplicates().set_index("sample_name")[['combinatorial_barcode', '384WP_well', 'Donor_ID']]
c.columns = ['combinatorial_barcode', 'plate_well', 'donor_id']
c.to_csv(
    os.path.join("metadata", "sciRNA-seq.PD192_PBMCs_196M.oligos_2019-09-20.csv"))

c = annot3.set_index('sample_name')[['combinatorial_barcode', '384WP_well', 'Cell_line']]
c.columns = ['combinatorial_barcode', 'plate_well', 'cell_line']
c.to_csv(
    os.path.join("metadata", "sciRNA-seq.PD193_fivelines_383k.oligos_2019-09-20.csv"))

c = annot4.set_index('sample_name')[['combinatorial_barcode', '384WP_well']]
c.columns = ['combinatorial_barcode', 'plate_well']
c.to_csv(
    os.path.join("metadata", "sciRNA-seq.PD193_humanmouse_765k.oligos_2019-09-20.csv"))

c = annot5.set_index('sample_name')[['combinatorial_barcode', '384WP_well']]
c.columns = ['combinatorial_barcode', 'plate_well']
c.to_csv(
    os.path.join("metadata", "sciRNA-seq.PD193_humanmouse_153M.oligos_2019-09-20.csv"))
