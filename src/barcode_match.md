Generate small file for testing
```bash
INPUT=$MUX
N=100000
OUTPUT=BSF_0674_HL7L5DMXX_1.head_${N}.bam

samtools view -H $INPUT > scifi_multiplexed.header
samtools view $INPUT  | head -n $N > scifi_multiplexed.head_${N}.body
cat scifi_multiplexed.header scifi_multiplexed.head_${N}.body \
| samtools view -b - > $OUTPUT
```
Identify where different barcodes are in the multiplexed file
```bash
MUX=/scratch/users/dbarreca/private/custom_demux/scRNA/BSF_0674_HL7L5DMXX/BSF_0674_HL7L5DMXX_1.bam
samtools view -h $MUX | head -n 8
# mux:      READ_STRUCTURE=21T8B16B78T
MUX=/scratch/users/dbarreca/private/custom_demux/scRNA/BSF_0674_HL7L5DMXX/BSF_0674_HL7L5DMXX_1_samples/BSF_0674_HL7L5DMXX_1#PD192_PBMCs_196M_77_A23_02.bam
samtools view -h $DEMUX | head -n 8
# demux:    READ_STRUCTURE=8M13B8B16M78T
```

```python
import os
import pandas as pd
import pysam
import numpy as np

# Read up barcode references
bcds = pd.read_csv(
    os.path.join("metadata", "PD191-193_multiplexing_barcodes.tsv"), sep='\t')
i7s = bcds['multiplexing_barcode'].unique().tolist()
r1s = bcds['combinatorial_barcode'].unique().tolist()
r2s = pd.read_csv(
    os.path.join("metadata", "737K-cratac-v1.reverse_complement.csv"))
r2s = r2s['original'].tolist()

# Read BAM file
n = 100000
bam_file = os.path.join(f"BSF_0674_HL7L5DMXX_1.head_{n}.bam")
bam = pysam.AlignmentFile(bam_file, check_sq=False)

bam_i7_file = os.path.join(f"BSF_0674_HL7L5DMXX_1.head_{n}.not_i7.bam")
bam_i7 = pysam.AlignmentFile(filename=bam_i7_file, mode="wb", header=bam.header)
bam_any_file = os.path.join(f"BSF_0674_HL7L5DMXX_1.head_{n}.not_any.bam")
bam_any = pysam.AlignmentFile(filename=bam_any_file, mode="wb", header=bam.header)

# n = 1000
t = int(n / 2)
res = np.zeros((t, 3))
write_any = False
write_i7 = False
for i, read in enumerate(bam):
    if read.is_read2:
        for v, f in [('write_any', bam_any), ('write_i7', bam_i7)]:
            if locals()[v]:
                read.is_paired = False
                read.is_read2 = False
                read.mate_is_unmapped = False
                f.write(read)
                locals()[v] = False
        continue
    if i == n:
        break
    # umi = read.seq[:8]
    r1 = read.seq[8:]
    tag = read.get_tag('BC')
    i7 = tag[:8]
    r2 = tag[8:]
    res[int(i / 2), :] = i7 in i7s, r1 in r1s, r2 in r2s
    if not (i7 in i7s):
        write_i7 = True
    if not res[int(i / 2), :].all():
        write_any = True
bam_i7.close()
bam_any.close()
```

Let's map the reads not matching to PhiX

```bash
# Get reference and build bowtie index
mkdir -p reference/phix/

wget -O reference/phix/phix.fa \
"https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1?report=fasta&log$=seqview&format=text"

bowtie2-build reference/phix/phix.fa reference/phix/phix

# Map with Bowtie2
for F in not_i7 not_any; do
bedtools bamtofastq \
-i BSF_0674_HL7L5DMXX_1.head_100000.${F}.bam \
-fq BSF_0674_HL7L5DMXX_1.head_100000.${F}.fastq \

bowtie2 \
-x reference/phix/phix \
-U BSF_0674_HL7L5DMXX_1.head_100000.${F}.fastq \
2> BSF_0674_HL7L5DMXX_1.head_100000.${F}.mapping.txt \
| samtools sort \
-o BSF_0674_HL7L5DMXX_1.head_100000.${F}.bowtie2.sorted.bam

rm BSF_0674_HL7L5DMXX_1.head_100000.${F}.fastq

done

```
Now let's update the estimates with the ammount of PhiX spike in in the samples
```python

c = res.all(1).sum()
v = [
    1 - (c / t),
    res.any(1).sum() / t,
    c / t, res[:, 0].sum() / t,
    res[:, 1].sum() / t,
    res[:, 2].sum() / t,
    1 - (res[:, [1, 2]].all(1).sum() / t),
    res[:, [1, 2]].any(1).sum() / t,
    res[:, [1, 2]].all(1).sum() / t]
k = [
    "None correct", "Any correct", "All correct",
    "i7 correct", "r1 correct", "r2 correct",
    "r1+r2, none correct", "r1+r2, any correct", "r1+r2, all correct"]
s = pd.Series(v, index=k)
# Assuming error rate of 0.25% per basepair(from PhiX spike-in)
s['All correct, seq error accounted'] = s['All correct'] + len(i7 + r1 + r2) * 0.0025
s['r1+r2, all correct, seq error accounted'] = s['r1+r2, all correct'] + len(i7 + r1 + r2) * 0.0025

# Account for PhIX
with open("BSF_0674_HL7L5DMXX_1.head_100000.not_any.mapping.txt", "r") as handle:
    content = handle.readlines()
s['phix_alignment'] = float(content[-1].split("%")[0]) / 100
s['phix_rate'] = (1 - s['r1+r2, all correct, seq error accounted']) * s['phix_alignment']
s['All correct, PhiX + seq error accounted'] = s['All correct, seq error accounted'] + s['phix_rate']
s['r1+r2, all correct, PhiX + seq error accounted'] = s['r1+r2, all correct, seq error accounted'] + s['phix_rate']

s.to_csv(f"BSF_0674_HL7L5DMXX_1.head_{n}.stats.csv", index=True, header=False)
```
