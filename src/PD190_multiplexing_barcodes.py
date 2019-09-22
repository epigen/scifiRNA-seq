import os
import pandas as pd

d = pd.read_excel(
    os.path.join(
        "metadata", "original", "PD000190_sixlines_round1_index_annotation.xlsx"
    )
)
d.loc[:, "combinatorial_barcode"] = (
    d["fixed_base1_N"] + d["index_seq"] + d["fixed_base2_V"]
)


# Add cell line mixing experiment
seqs = [("GGACTCCT", "01"), ("TAGGCATG", "02"), ("CTCTCTAC", "03"), ("CGAGGCTG", "04")]
annot = pd.concat(
    [d.assign(multiplexing_barcode=b, multiplexing_number=i) for b, i in seqs]
).reset_index(drop=True)
annot = annot.assign(
    sample_name="PD190_sixlines"
    + "_"
    + annot["cell_line"].str.replace("-", "")
    + "_"
    + annot["plate_well"]
    + "_"
    + annot["multiplexing_number"]
)

# Add species mixing experiment
seqs = [("TAAGGCGA", "01"), ("CGTACTAG", "02"), ("AGGCAGAA", "03"), ("TCCTGAGC", "04")]
annot2 = pd.concat(
    [d.assign(multiplexing_barcode=b, multiplexing_number=i) for b, i in seqs]
).reset_index(drop=True)
annot2 = annot2.assign(
    sample_name="PD190_humanmouse"
    + "_"
    + annot["plate_well"]
    + "_"
    + annot["multiplexing_number"]
)

annot = annot.append(annot2)
annot[["sample_name", "multiplexing_barcode", "combinatorial_barcode"]].to_csv(
    os.path.join("metadata", "PD190_multiplexing_barcodes.tsv"), sep="\t", index=False
)

annot.set_index("sample_name").to_csv(
    os.path.join("metadata", "sciRNA-seq.PD190.oligos_2019-09-05.full_annotation.csv")
)


c = annot.loc[annot['sample_name'].str.contains("humanmouse"), :]
c = c.drop(["multiplexing_number", "multiplexing_barcode"], axis=1)
c['sample_name'] = c['sample_name'].str.replace(r"_\d\d$", "")
c.drop_duplicates().set_index("sample_name").to_csv(
    os.path.join("metadata", "sciRNA-seq.PD190_humanmouse.oligos_2019-09-05.csv")
)

c = annot.loc[annot['sample_name'].str.contains("sixlines"), :]
c = c.drop(["multiplexing_number", "multiplexing_barcode"], axis=1)
c['sample_name'] = c['sample_name'].str.replace(r"_\d\d$", "")
c.drop_duplicates().set_index("sample_name").to_csv(
    os.path.join("metadata", "sciRNA-seq.PD190_sixlines.oligos_2019-09-05.csv")
)
