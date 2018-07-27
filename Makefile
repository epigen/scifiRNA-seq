.DEFAULT_GOAL := run

run:
	sh src/process_barcodes.sh
	sh src/inspect_barcodes.py

all: run

.PHONY: run all
