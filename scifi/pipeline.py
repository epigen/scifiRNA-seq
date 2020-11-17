#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
The entry point to running the scifi_pipeline on several samples.
"""

import os
import argparse

import pandas as pd

from scifi.map import map_command
from scifi.filter import filter_command
from scifi.join import join_command
from scifi import _LOGGER, setup_config


SAMPLE_ANNOTATION_COLUMNS = [
    "sample_name",
    "annotation",
    "variables",
    "species_mixing",
    "expected_cell_number",
]


def build_cli() -> argparse.ArgumentParser:
    _LOGGER.debug("Setting up CLI parser.")
    parser = argparse.ArgumentParser()

    sp = parser.add_subparsers(dest="command", required=True)
    all_cmd = sp.add_parser("all")
    map_cmd = sp.add_parser("map")
    filter_cmd = sp.add_parser("tracks")
    filter_cmd = sp.add_parser("filter")
    join_cmd = sp.add_parser("join")
    report_cmd = sp.add_parser("report")
    # Arguments common to all commands
    for cmd in [all_cmd, map_cmd, filter_cmd, join_cmd, report_cmd]:
        if cmd != join_cmd:
            _help = "Whether to use job arrays (only supported for SLURM)."
            cmd.add_argument("--arrayed", help=_help, action="store_true")
            _help = "Size of job arrays."
            cmd.add_argument("--array-size", dest="array_size", help=_help, type=int)
        _help = "Whether to not submit any jobs."
        cmd.add_argument("-d", "--dry-run", action="store_true", help=_help)
        _help = "Whether to only run samples marked with '1' in toggle column."
        cmd.add_argument("-t", "--toggle", action="store_true", help=_help)
        _help = "Samples to subset. Comma delimited."
        cmd.add_argument(
            "-s", "--sample-subset", dest="sample_subset", help=_help, default="",
        )
        _help = (
            "YAML configuration file."
            " If not provided will use the one distributed"
            " with the package, or one in a file in"
            " ~/.scifiRNA-seq.config.yaml"
        )
        cmd.add_argument("-c", "--config-file", dest="config_file", help=_help, default=None)
        _help = (
            "Directory to output files."
            "By default it will be 'pipeline_output'"
            " in current directory"
        )
        _def = os.path.join(os.path.curdir, "pipeline_output")
        cmd.add_argument(
            "-o", "--output-dir", dest="root_output_dir", help=_help, default=_def,
        )
        _help = (
            "CSV file with sample annotation. One row per sample."
            " Mandatory columns are: " + ", ".join(SAMPLE_ANNOTATION_COLUMNS) + "."
        )
        cmd.add_argument(dest="sample_annotation", help=_help)
        # Requirements for a sample (rows of `sample_annotation`):
        # - sample_name: str (a name for a sample)
        # - annotation: str (a CSV file with mapping of round1 barcodes and their attributes)
        # # In its turn, `annotation` should have the following columns:
        # # - sample_name: str (name of round1)
        # # - combinatorial_barcode: str (DNA sequence of round1 barcode)
        # # - plate_well: str (name of round1 well)
        # # - ...: str, optional (additional attributes, should be listed in `variables` below)
        # - variables: str (which of the columns of `annotation` should be used to annotate round1 wells)
        # - species_mixing: bool (whether the experiment is a species mixing)
        # - expected_cell_number: int (expected cell number)
        # - toggle: bool, optional (whether the sample should be run)
        # - pass_qc: bool, optional (whether the sample passed QC and should be run)

    # Map-specific
    _help = (
        "A path or glob to input BAM files."
        " This can include variables in curly brackets that will be replaced "
        " with values from the sample metadata "
        "(e.g. {sample_name}) or round1 metadata (e.g. {plate_well})."
        " Example: /lab/seq/{flowcell}/{flowcell}#*_{sample_name}.bam"
    )
    map_cmd.add_argument("--input-bam-glob", dest="input_bam_glob", help=_help, required=True)

    # Filter-specific
    _help = "Whether to correct round2 barcodes."
    filter_cmd.add_argument(
        "--correct-r2-barcodes", dest="correct_r2_barcodes", action="store_true", help=_help,
    )
    _help = "Whitelist of round2 barcodes."
    filter_cmd.add_argument("--r2-barcodes-whitelist", dest="r2_barcodes_whitelist", help=_help)
    _help = "Whether to overwrite exitsing files."
    filter_cmd.add_argument("--overwrite", dest="overwrite", action="store_true", help=_help)
    return parser


def main(cli=None):
    _LOGGER.info("scifi-RNA-seq pipeline")

    # Parse arguments and config
    args = build_cli().parse_args(cli)
    args.sample_subset = args.sample_subset.split(",") if args.sample_subset != "" else []
    _LOGGER.debug(args)
    _CONFIG = setup_config(args.config_file)
    _LOGGER.debug(_CONFIG)

    # Read and prepare sample annotation sheet
    df = pd.read_csv(args.sample_annotation)
    if args.toggle and "toggle" in df.columns:
        df = df.query("toggle == 1")
    if args.sample_subset:
        df = df.loc[df["sample_name"].isin(args.sample_subset), :]

    if df.empty:
        _LOGGER.error("No samples have been selected. Terminating.")
        return 1
    # # assume that species mixing is False is missing
    df.loc[:, "species_mixing"] = df["species_mixing"].fillna(value=0).astype(bool)
    # # leave expected cell numbers to default in respective functions
    df.loc[:, "expected_cell_number"] = df["expected_cell_number"].fillna(value=200000)

    s = "\n\t - " + "\n\t - ".join(df["sample_name"])
    _LOGGER.info(f"Samples to submit:{s}")

    # Process per sample
    for i, sample in df.iterrows():
        _LOGGER.debug(f"Doing sample {sample['sample_name']}")

        sample_name = sample["sample_name"]
        r1_annotation_file = sample["annotation"]
        r1_annotation = pd.read_csv(r1_annotation_file).set_index("sample_name")
        sample_out_dir = os.path.join(args.root_output_dir, sample_name)
        os.makedirs(sample_out_dir, exist_ok=True)
        sample_attributes = sample["variables"].split(",")

        if args.command in ["all", "map"]:
            _LOGGER.debug(f"Running map command with sample {sample_name}")
            map_command(args, sample_name, sample_out_dir, r1_annotation, _CONFIG)
        if args.command in ["all", "tracks"]:
            _LOGGER.warning("The tracks step is not yet implemented. Skipping.")
            pass
            # _LOGGER.debug(f"Running tracks command with sample {sample_name}")
            # tracks_command(args)
        if args.command in ["all", "filter"]:
            _LOGGER.debug(f"Running filter command with sample {sample_name}")
            filter_command(
                args,
                _CONFIG,
                sample_name,
                sample_out_dir,
                r1_annotation,
                r1_annotation_file=r1_annotation_file,
                r1_attributes=sample_attributes,
                species_mixture=sample["species_mixing"],
                expected_cell_number=sample["expected_cell_number"],
                correct_r2_barcodes=False,
                correct_r2_barcodes_file=None,
                dry_run=args.dry_run,
            )
        if args.command in ["all", "join"]:
            _LOGGER.debug(f"Running join command with sample {sample_name}")
            join_command(
                sample_name,
                sample_out_dir,
                r1_attributes=sample_attributes,
                species_mixture=sample["species_mixing"],
                correct_r2_barcodes=False,
            )
        if args.command in ["all", "report"]:
            _LOGGER.warning("The report step is not yet implemented. Skipping.")
            pass
            # _LOGGER.debug(f"Running report command with sample {sample_name}")
            # report_command(args)
