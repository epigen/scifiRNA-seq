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
EXPECTED_CELL_NUMBER ?= 200000
MIN_UMI_OUTPUT ?= 3
VARIABLES ?= "plate_well"
ARRAY_SIZE ?= 24
SPECIES_MIXING ?= 1
SPECIES_MIX_FLAG := 
ifeq ($(SPECIES_MIXING), 1)
	SPECIES_MIX_FLAG := --species-mixture
	STAR_DIR := /home/arendeiro/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/indexed_STAR_2.7.0e/
	GTF_FILE := /home/arendeiro/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/Homo_sapiens-Mus_musculus.Ensembl92.dna.primary_assembly.Tcr_lambda_spiked.gtf
endif
CHUNKS ?= 1000
CHUNK_BATCH_SIZE ?= 25


trim: parse
	# only for gRNA samples
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
	--cpus=4 \
	--mem=60000 \
	--queue=shortq \
	--time=08:00:00 \
	--array-size=$(ARRAY_SIZE) \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--star-exe=$(STAR_EXE) \
	--star-dir=$(STAR_DIR) \
	--gtf=$(GTF_FILE)
	$(info "scifi_pipeline: done")

maketracks: parse
	$(info "scifi_pipeline: maketracks")
	sh src/scifi_pipeline.maketracks.sh \
	--run-name=$(RUN_NAME) \
	--annotation=$(ANNOTATION) \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--cpus=1 \
	--mem=8000 \
	--queue=shortq \
	--time=01:00:00 \
	--array-size=$(ARRAY_SIZE)
	$(info "scifi_pipeline: done")

filter: parse
	$(info "scifi_pipeline: filter")
	sh src/scifi_pipeline.filter.sh \
	--run-name=$(RUN_NAME) \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--annotation=$(ANNOTATION) \
	--expected-cell-number=$(EXPECTED_CELL_NUMBER) \
	--min-umi-output=${MIN_UMI_OUTPUT} \
	--variables=$(VARIABLES) \
	--species-mixture=$(SPECIES_MIXING) \
	--cpus=1 \
	--mem=48000 \
	--queue=shortq \
	--time=06:00:00 \
	--array-size=$(ARRAY_SIZE)
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

barcodes: parse
	$(info "scifi_pipeline: barcodes")
	for CHUNK in $(shell seq 0 $(CHUNK_BATCH_SIZE) $(CHUNKS)); \
	do \
		sbatch -J scifi_pipeline.barcodes.$(RUN_NAME).$(CHUNKS).$${CHUNK} \
		-o $(ROOT_OUTPUT_DIR)/scifi_pipeline.barcodes.$(RUN_NAME).$(CHUNKS).$${CHUNK}.log \
		-p shortq --mem 16000 --cpus 1 --time 04:00:00 \
		--wrap "python3 -u src/scifi_pipeline.barcode_matcher.py \
		--metrics $(ROOT_OUTPUT_DIR)/$(RUN_NAME).metrics.csv.gz \
		--output-prefix $(ROOT_OUTPUT_DIR)/$(RUN_NAME). \
		--total-chunk $(CHUNKS) \
		--start-chunk $${CHUNK} \
		--batch-size $(CHUNK_BATCH_SIZE)"; \
	done

fix: parse
	cat $(ROOT_OUTPUT_DIR)/$(RUN_NAME).fixed_barcodes.*-*.tsv | \
	sort > $(ROOT_OUTPUT_DIR)/$(RUN_NAME).fixed_barcodes.mapping.tsv

filter_corrected: parse
	$(info "scifi_pipeline: filter_corrected")
	sh src/scifi_pipeline.filter.sh \
	--run-name=$(RUN_NAME) \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--annotation=$(ANNOTATION) \
	--expected-cell-number=$(EXPECTED_CELL_NUMBER) \
	--min-umi-output=${MIN_UMI_OUTPUT} \
	--variables=$(VARIABLES) \
	--species-mixture=$(SPECIES_MIXING) \
	--cpus=1 \
	--mem=8000 \
	--queue=shortq \
	--time=01:00:00 \
	--correct-barcodes=1 \
	--correct-barcode-file=$(ROOT_OUTPUT_DIR)/$(RUN_NAME).fixed_barcodes.mapping.tsv \
	--array-size=$(ARRAY_SIZE)
	$(info "scifi_pipeline: done")

join_corrected: parse
	$(info "scifi_pipeline: join_corrected")
	sh src/scifi_pipeline.join.sh \
	--run-name=$(RUN_NAME) \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--variables=$(VARIABLES) \
	--species-mixture=$(SPECIES_MIXING) \
	--cpus=1 \
	--mem=12000 \
	--queue=shortq \
	--time=08:00:00 \
	--correct-barcodes=1
	$(info "scifi_pipeline: done")

report: parse
	$(info "scifi_pipeline: report")

	mkdir -p results/$(RUN_NAME)

	sbatch -J scifi_pipeline.report.$(RUN_NAME) \
	-o $(ROOT_OUTPUT_DIR)/scifi_pipeline.report.log \
	-p longq --mem 120000 --cpus 4 --time 3-00:00:00 \
	--wrap "python3 -u src/scifi_pipeline.report.py \
	$(ROOT_OUTPUT_DIR)/$(RUN_NAME).metrics.csv.gz \
	results/$(RUN_NAME)/$(RUN_NAME). \
	--plotting-attributes $(VARIABLES) $(SPECIES_MIX_FLAG)"

	sbatch -J scifi_pipeline.report-exon.$(RUN_NAME) \
	-o $(ROOT_OUTPUT_DIR)/scifi_pipeline.report-exon.log \
	-p longq --mem 120000 --cpus 4 --time 3-00:00:00 \
	--wrap "python3 -u src/scifi_pipeline.report.py \
	$(ROOT_OUTPUT_DIR)/$(RUN_NAME).exon.metrics.csv.gz \
	results/$(RUN_NAME)/$(RUN_NAME).exon. \
	--plotting-attributes $(VARIABLES) $(SPECIES_MIX_FLAG)"
	$(info "scifi_pipeline: done")

	sbatch -J scifi_pipeline.report_corrected.$(RUN_NAME) \
	-o $(ROOT_OUTPUT_DIR)/scifi_pipeline.report_corrected.log \
	-p longq --mem 120000 --cpus 4 --time 3-00:00:00 \
	--wrap "python3 -u src/scifi_pipeline.report.py \
	$(ROOT_OUTPUT_DIR)/$(RUN_NAME).metrics_corrected.csv.gz \
	results/$(RUN_NAME)/$(RUN_NAME).corrected. \
	--plotting-attributes $(VARIABLES) $(SPECIES_MIX_FLAG)"

	sbatch -J scifi_pipeline.report_corrected-exon.$(RUN_NAME) \
	-o $(ROOT_OUTPUT_DIR)/scifi_pipeline.report_corrected-exon.log \
	-p longq --mem 120000 --cpus 4 --time 3-00:00:00 \
	--wrap "python3 -u src/scifi_pipeline.report.py \
	$(ROOT_OUTPUT_DIR)/$(RUN_NAME).exon.metrics_corrected.csv.gz \
	results/$(RUN_NAME)/$(RUN_NAME).exon.corrected. \
	--plotting-attributes $(VARIABLES) $(SPECIES_MIX_FLAG)"
	$(info "scifi_pipeline: done")


matches:
	$(info "scifi: barcode_match")
	for N in 1000 5000 10000 20000 200000 2000000 20000000 200000000; \
	do \
		sbatch -J scifi.barcode_match.$${N} \
		-o results/barcode_matches/scifi.barcode_match.$${N}.log \
		-p mediumq --mem 32000 --cpus 1 --time 1-08:00:00 \
		--wrap "python3 -u src/barcode_profiling.py \
		-n $${N} metadata/annotation.csv"; \
	done


metrics: matches



all: trim map filter join barcodes fix filter_corrected join_corrected report

clean:
	find . -name "*bam" ! -name external -delete

.PHONY: trim map filter join barcodes fix filter_corrected join_corrected report clean
.SILENT: all
