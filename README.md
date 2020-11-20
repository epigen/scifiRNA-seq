> :warning: **This repository is experimental**: Testing is very limited and we
provide no promise of continually developing or maintain it!

# scifiRNA-seq pipeline

This is a pipeline to process
[scifi-RNA-seq](https://www.biorxiv.org/content/10.1101/2019.12.17.879304v1)
data.

This repository was **not** used to process the data in the publication (source
code will be made available upon publication), but is meant to be a more
portable and performant version of that.

## Requirements

This pipeline has been developed only with a Linux system in mind and other
are not supported, feel free to contribute a PR if facing problems in other
systems.

- Python >=3.7
- SLURM (not a hard requirement but best supported currently)

The Python requirements can be seen in the [requirements.txt](requirements.txt)
file and will be installed automatically by ``pip``.

### Additional requirements:

 - samtools
 - [STAR aligner](https://github.com/alexdobin/STAR) (we used version 2.7.0e and
recommend using that or above).
 - [featureCounts](http://subread.sourceforge.net/)


This pipeline was made to take advantage of high parallelization with a high
performance computing cluster. It currently only supports the SLURM job
scheduler and uses job arrays.

## Installation
In due time the pipeline will be uploaded to Pypi, but for now either set up git
with SSH and use (requires repository access):

```bash
pip install git+ssh://git@github.com/epigen/scifiRNA-seq.git
```
Or clone the repository and install it:

```bash
git clone https://github.com/epigen/scifiRNA-seq
cd scifiRNA-seq
pip install -e .
```

Make sure to use an up-to-date ``pip`` version.

After installation, an executable ``scifi`` will be added to the user's local
bin directory (usually `~/.local/bin`), make sure that is in your `PATH`
variable.

## Configuration and logging

scifi pipeline ships with a
[default configuration file](scifi/config/default.yaml).


In the configuration file, location of software dependencies (STAR, featureCounts), but also location of static files (e.g. genome index, GTF file) can be specified.

To configure the pipeline to a specific environment, write a file with the same structure to `~/.scifi.config.yaml` to avoid passing a specific configuration at runtime repeatedly (but possible with the -c option).

A log file is written to `~/.scifi.log.txt` in addition to the command-line
output.


## Running the pipeline

### Required input files

#### BAM files

The pipeline expects unaligned BAM input. These files should be produced by
demultiplexing the raw base calls based on any sample index and the round1
index, yielding one file per experiment, per round1 barcode.

Each read should have a the following tags:
 - `BC`: a concatenation of the sample barcode and the round1 barcode (22 bp);
 - `RX`: the UMI barcode (8bp)

We use a [custom fork of Picard tools](https://github.com/DanieleBarreca/picard)
for demultiplexing. The commands used for demultiplexing will be made available
in the repository of the publication.

#### Sample annotation

To run the pipeline prepare a CSV annotation of your samples.
Each row represents one experiment.

Mandatory columns are:
 - `sample_name`: A name for the sample
 - `annotation`: CSV file with annotation for each round1 well.
 - `variables`:  Variables in CSV `annotation` file above to use. A
comma-delimited string with various values (no spaces). Use software capable of appropriately quoting fields to produce CSV file.
 - `species_mixing`: Whether the experiment is a Banyard experiment with cells
from two organisms. Use `1` or `0` for this column.
 - `expected_cell_number`: The expected number of cells. This has no influence
on the actual number of reported cells but is used as a comparison to that.

A special boolean column `toggle` can be used to select/deselect samples at runtime.
Any other columns may be added but will not be used by the pipeline.

Example:

```csv
sample_name,toggle,protocol,batch,cell_line,expected_cell_number,material,organism,species_mixing,flowcell,variables
scifi-RNA-seq_PDXYZ_SAMPLE_X_200k-nuclei,0,scifiRNA-seq,SCI004,Jurkat,200000,nuclei,human,1,BSF_XXXXX,"plate_well,treatment,knockout"
scifi-RNA-seq_PDXYZ_SAMPLE_X_1M-nuclei,1,scifiRNA-seq,SCI004,Jurkat,1000000,nuclei,human,1,BSF_XXXXX,"plate_well,treatment,knockout"
```

#### Round1 plate well annotation

A CSV file with one row per well.

Mandatory columns are:
 - `sample_name`: A name for this combination of sample and well;
 - `plate_well`: The plate well code e.g. A01, F08;
 - `combinatorial_barcode`: the sequence of the round1 barcode;

Supplement this file with any additional columns to for example annotate
experimental conditions. Add the name of those columns to the `variables` field
of the sample CSV annotation file in order to have them used by the pipeline.

Example:

```csv
sample_name,combinatorial_barcode,plate_well,treatment,knockout
scifi-RNA-seq_PDXYZ_SAMPLE_X_1M-nuclei_A01,AAGTGATTAGCAA,A01,DMSO,knockoutA
scifi-RNA-seq_PDXYZ_SAMPLE_X_1M-nuclei_A03,AGAATCCCCCTAA,A03,DMSO,knockoutB
scifi-RNA-seq_PDXYZ_SAMPLE_X_1M-nuclei_A05,ACCTGGGAAACTA,A05,Gefitinib,knockoutA
scifi-RNA-seq_PDXYZ_SAMPLE_X_1M-nuclei_A07,ATACCTCCCAGGA,A07,Gefitinib,knockoutB
```

### Runnning

The ``scifi`` executable is placed in the path `pip` installs software to. In linux systems this will usually be `~/.local/bin`.
In order to call the command without refering to that location, add `~/.local/bin` to you `$PATH` variable.

To see the help for the pipeline:
```bash
scifi --help
```

The pipeline has several commands. To see the help for a specific command:
```bash
scifi map --help
```

To run a command for all samples simply run:
```bash
scifi \
    map \
        --input-bam-glob /lab/seq/{flowcell}/{flowcell}#*_{sample_name}.bam \
        metadata/annotation.csv
```

A new configuration file can be passed at runtime with the "-c" option or
values specfied in `~/.scifi.config.yaml`.

If not using SLURM, provide a value in the configuration unde `submission_command` that is the command to be called to execute the job e.g. `sh`.

A dry run is possible with the option `-d/--dry-run`, which will produce job files to be executed - useful for debugging or running in a particular way of choice.


## Pipeline outputs

The most relevant outputs include:
 - Per well, mapped and gene tagged BAM file;
 - CSV file with summary statistics per barcode;
 - CSV file with expression values per cell, per gene;
 - h5ad gene expression file;

Additional outputs include various visualizations related to graphics presented
in the preprint, such as "knee plots" and species mixture plots.

The pipeline is a little wasteful in that it trades disk space usage for speed.
If space is a limiting factor we recommend deleting aligned BAM files after a
successful run.


## Contributions

- Andre F. Rendeiro
