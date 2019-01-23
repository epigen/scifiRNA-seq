.DEFAULT_GOAL := run

run:
	sh src/main_pipeline.sh  # this actualy should be run manually
	python -u src/plot_performance.SCI016.py

all: run

.PHONY: run all
