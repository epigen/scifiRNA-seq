.DEFAULT_GOAL := run

run:
	sh src/process_barcodes.sh
	python -u src/inspect_barcodes.py

all: run

.PHONY: run all
