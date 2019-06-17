.DEFAULT_GOAL := all

# Change these parameters
RUN_NAME       := SCI024_Tcell_s
FLOWCELL       := BSF_0624_HJY5CDMXX
N_LANES        := 2
ANNOTATION     := /scratch/lab_bock/shared/projects/sci-rna/metadata/sciRNA-seq.SCI024.oligos_2019-05-17.csv
ROOT_OUTPUT_DIR:= /scratch/lab_bock/shared/projects/sci-rna/data/$(RUN_NAME)
STAR_EXE       := /home/arendeiro/workspace/STAR-2.7.0e/bin/Linux_x86_64_static/STAR
STAR_DIR       := /home/arendeiro/resources/genomes/hg38/indexed_STAR-2.7.0e/
GTF_FILE       := /home/arendeiro/resources/genomes/hg38/10X/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf

# Don't change below unless you need different resources
map:
	echo "scifi_pipeline: map"
	sh src/scifi_pipeline.map.sh \
	--run-name=$(RUN_NAME) \
	--flowcell=$(FLOWCELL) \
	--n-lanes=$(N_LANES) \
	--annotation=$(ANNOTATION) \
	--cpus=4 \
	--mem=50000 \
	--queue=shortq \
	--time=08:00:00 \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--star-exe=$(STAR_EXE) \
	--star-dir=$(STAR_DIR) \
	--gtf=$(GTF_FILE)

filter:
	echo "scifi_pipeline: filter"
	sh src/scifi_pipeline.filter.sh \
	--run-name=$(RUN_NAME) \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--cpus=1 \
	--mem=8000 \
	--queue=shortq \
	--time=08:00:00

join:
	echo "scifi_pipeline: join"
	sh src/scifi_pipeline.join.sh \
	--run-name=$(RUN_NAME) \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--cpus=1 \
	--mem=12000 \
	--queue=shortq \
	--time=08:00:00

plot:
	echo "scifi_pipeline: plot"
	python3 -u src/sci-RNA-seq.report.py \
	$(ROOT_OUTPUT_DIR)/$(RUN_NAME).metrics.csv.gz \
	results/$(RUN_NAME). \
	--plotting-attributes plate donor_id

	python3 -u src/sci-RNA-seq.report.py \
	$(ROOT_OUTPUT_DIR)/$(RUN_NAME).exon.metrics.csv.gz \
	results/$(RUN_NAME).exon. \
	--plotting-attributes plate donor_id

all: map filter join plot

clean:
	find . -name "*bam" -delete

.PHONY: map filter join plot clean
.SILENT: all
