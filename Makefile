.DEFAULT_GOAL := all

# Change these parameters
STAR_EXE       := /home/arendeiro/workspace/STAR-2.7.0e/bin/Linux_x86_64_static/STAR
STAR_DIR       := /home/arendeiro/resources/genomes/hg38/indexed_STAR-2.7.0e/
STAR_DIR       := /home/arendeiro/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/indexed_STAR_2.7.0e/
GTF_FILE       := /home/arendeiro/resources/genomes/hg38/10X/refdata-cellranger-GRCh38-1.2.0/genes/genes.gtf
GTF_FILE       := /home/arendeiro/resources/genomes/hg38_mm10_transgenes_Tcrlibrary/Homo_sapiens-Mus_musculus.Ensembl92.dna.primary_assembly.Tcr_lambda_spiked.gtf


trim: parse
	@echo "scifi_pipeline: trim"
	sh src/scifi_pipeline.trim.sh \
	--run-name=$(RUN_NAME) \
	--flowcell=$(FLOWCELL) \
	--n-lanes=$(N_LANES) \
	--annotation=$(ANNOTATION) \
	--cpus=1 \
	--mem=8000 \
	--queue=shortq \
	--time=08:00:00 \
	--output-dir=$(ROOT_OUTPUT_DIR)

map: parse
	@echo "scifi_pipeline: map"
	sh src/scifi_pipeline.map.sh \
	--run-name=$(RUN_NAME) \
	--flowcell=$(FLOWCELL) \
	--n-lanes=$(N_LANES) \
	--annotation=$(ANNOTATION) \
	--cpus=4 \
	--mem=80000 \
	--queue=shortq \
	--time=08:00:00 \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--star-exe=$(STAR_EXE) \
	--star-dir=$(STAR_DIR) \
	--gtf=$(GTF_FILE)

filter: parse
	@echo "scifi_pipeline: filter"
	sh src/scifi_pipeline.filter.sh \
	--run-name=$(RUN_NAME) \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--annotation=$(ANNOTATION) \
	--cpus=1 \
	--mem=8000 \
	--queue=shortq \
	--time=08:00:00

join: parse
	@echo "scifi_pipeline: join"
	sh src/scifi_pipeline.join.sh \
	--run-name=$(RUN_NAME) \
	--output-dir=$(ROOT_OUTPUT_DIR) \
	--cpus=1 \
	--mem=12000 \
	--queue=shortq \
	--time=08:00:00

plot: parse
	@echo "scifi_pipeline: plot"
	sbatch -J scifiRNA-seq.$(RUN_NAME).plot \
	-o $(ROOT_OUTPUT_DIR)/$(RUN_NAME).plot.log \
	-p shortq --mem 120000 --cpus 2 \
	--wrap "python3 -u src/scifi_pipeline.report.py \
	$(ROOT_OUTPUT_DIR)/$(RUN_NAME).metrics.csv.gz \
	results/$(RUN_NAME). \
	--plotting-attributes plate donor_id"

	sbatch -J scifiRNA-seq.$(RUN_NAME).plot-exon \
	-o $(ROOT_OUTPUT_DIR)/$(RUN_NAME).plot-exon.log \
	-p shortq --mem 120000 --cpus 2 \
	--wrap "python3 -u src/scifi_pipeline.report.py \
	$(ROOT_OUTPUT_DIR)/$(RUN_NAME).exon.metrics.csv.gz \
	results/$(RUN_NAME).exon. \
	--plotting-attributes plate donor_id"

parse:
	@[ "${RUN_NAME}" ] || ( echo "'RUN_NAME' is not set"; exit 1 )
	@[ "${FLOWCELL}" ] || ( echo "'FLOWCELL' is not set"; exit 1 )
	@[ "${N_LANES}" ] || ( echo "'N_LANES' is not set"; exit 1 )
ifeq ($(origin ANNOTATION),undefined)
ANNOTATION := "/scratch/lab_bock/shared/projects/sci-rna/metadata/sciRNA-seq.SCI024.oligos_2019-05-17.csv"
endif
ifeq ($(origin ROOT_OUTPUT_DIR),undefined)
ROOT_OUTPUT_DIR := /scratch/lab_bock/shared/projects/sci-rna/data/$(RUN_NAME)
endif


all: trim map filter join plot

clean:
	find . -name "*bam" -delete

.PHONY: trim map filter join plot clean
.SILENT: all
