.DEFAULT_GOAL := all

# These parameters can be overwritten
STAR_EXE := /home/arendeiro/workspace/STAR-2.7.0e/bin/Linux_x86_64_static/STAR
STAR_DIR := /home/arendeiro/resources/genomes/hg38/indexed_STAR-2.7.0e/
GTF_FILE := /home/arendeiro/resources/genomes/hg38/10X/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf

parse:
	@[ "${RUN_NAME}" ] || ( echo "'RUN_NAME' is not set"; exit 1 )
	@[ "${FLOWCELL}" ] || ( echo "'FLOWCELL' is not set"; exit 1 )
	@[ "${N_LANES}" ] || ( echo "'N_LANES' is not set"; exit 1 )
	@[ "${N_BARCODES}" ] || ( echo "'N_BARCODES' is not set"; exit 1 )
ANNOTATION ?= "/scratch/lab_bock/shared/projects/sci-rna/metadata/sciRNA-seq.PD190_humanmouse.oligos_2019-09-05.csv"
ROOT_OUTPUT_DIR ?= /scratch/lab_bock/shared/projects/sci-rna/data/$(RUN_NAME)
VARIABLES ?= "plate_well"
SPECIES_MIXING ?= 1
SP := 
ifeq ($(SPECIES_MIXING), 1)
	SP := --species-mixture
	STAR_DIR := /home/arendeiro/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/indexed_STAR_2.7.0e/
	GTF_FILE := /home/arendeiro/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/Homo_sapiens-Mus_musculus.Ensembl92.dna.primary_assembly.Tcr_lambda_spiked.gtf
endif

trim: parse
	$(info "scifi_pipeline: trim")
	sh src/scifi_pipeline.trim.sh \
	--run-name=$(RUN_NAME) \
	--flowcell=$(FLOWCELL) \
	--n-lanes=$(N_LANES) \
	--annotation=$(ANNOTATION) \
	--cpus=1 \
	--mem=8000 \
	--queue=longq \
	--time=08:00:00 \
	--output-dir=$(ROOT_OUTPUT_DIR)
	$(info "scifi_pipeline: done")

map: parse
	$(info "scifi_pipeline: map")
	sh src/scifi_pipeline.map.sh \
	--run-name=$(RUN_NAME) \
	--flowcell=$(FLOWCELL) \
	--n-lanes=$(N_LANES) \
	--n-barcodes=$(N_BARCODES) \
	--annotation=$(ANNOTATION) \
	--cpus=1 \
	--mem=80000 \
	--queue=longq \
	--time=00:30:00 \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--star-exe=$(STAR_EXE) \
	--star-dir=$(STAR_DIR) \
	--gtf=$(GTF_FILE)
	$(info "scifi_pipeline: done")

filter: parse
	$(info "scifi_pipeline: filter")
	sh src/scifi_pipeline.filter.sh \
	--run-name=$(RUN_NAME) \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--annotation=$(ANNOTATION) \
	--variables=$(VARIABLES) \
	--species-mixture=$(SPECIES_MIXING) \
	--cpus=1 \
	--mem=8000 \
	--queue=longq \
	--time=08:00:00
	$(info "scifi_pipeline: done")

join: parse
	$(info "scifi_pipeline: join")
	sh src/scifi_pipeline.join.sh \
	--run-name=$(RUN_NAME) \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--variables=$(VARIABLES) \
	--species-mixture=$(SPECIES_MIXING) \
	--cpus=1 \
	--mem=12000 \
	--queue=shortq \
	--time=08:00:00
	$(info "scifi_pipeline: done")

report: parse
	$(info "scifi_pipeline: report")

	sbatch -J scifiRNA-seq.$(RUN_NAME).report \
	-o $(ROOT_OUTPUT_DIR)/$(RUN_NAME).report.log \
	-p longq --mem 120000 --cpus 2 \
	--wrap "python3 -u src/scifi_pipeline.report.py \
	$(ROOT_OUTPUT_DIR)/$(RUN_NAME).metrics.csv.gz \
	results/$(RUN_NAME). \
	--plotting-attributes $(VARIABLES) $(SP)"

	sbatch -J scifiRNA-seq.$(RUN_NAME).report-exon \
	-o $(ROOT_OUTPUT_DIR)/$(RUN_NAME).report-exon.log \
	-p longq --mem 120000 --cpus 2 \
	--wrap "python3 -u src/scifi_pipeline.report.py \
	$(ROOT_OUTPUT_DIR)/$(RUN_NAME).exon.metrics.csv.gz \
	results/$(RUN_NAME).exon. \
	--plotting-attributes $(VARIABLES) $(SP)"
	$(info "scifi_pipeline: done")

all: trim map filter join report

clean:
	find . -name "*bam" -delete

.PHONY: trim map filter join report clean
.SILENT: all
