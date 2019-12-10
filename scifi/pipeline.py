#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
The entry point to running the scifi_pipeline on several samples.
"""

import os
from argparse import ArgumentParser

import yaml
import pandas as pd

from scifi.utils import ct
from scifi.map import map_command
from scifi.filter import filter_command


def build_cli():
    parser = ArgumentParser()

    sp = parser.add_subparsers(dest="command")
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
        cmd.add_argument("-s", "--sample-subset", help=_help, default="")
        _help = "Samples to subset. Comma delimited."
        cmd.add_argument(
            "-c", "--config-file", dest="config_file", help=_help, default="")
        _help = (
            "Directory to output files. By default it will be 'pipeline_output'"
            " in current directory")
        cmd.add_argument("-o", "--output-dir", dest="root_output_dir", help=_help)
        _help = (
            "YAML configuration file."
            " If not provided will use the one distributed"
            " with the package, or one in a file in"
            " ~/.scifiRNA-seq.config.yaml")
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

    _help = (
        "A path or glob to input BAM files."
        " This can include variables in curly brackets that will be replaced "
        " with values from "
        " the sample metadata (e.g. {sample_name}) or round1 metadata"
        " (e.g. {plate_well})."
        " Example: /scratch/lab/sequences/{flowcell}/{flowcell}#*_{sample_name}.bam")
    map_cmd.add_argument("--input-bam-glob", dest="input_bam_glob", help=_help)

    _help = "Whether to correct round2 barcodes."
    filter_cmd.add_argument("--correct-r2-barcodes", dest="correct_r2_barcodes", action="store_true", help=_help)
    _help = "Whitelist of round2 barcodes."
    filter_cmd.add_argument("--r2-barcodes-whitelist", dest="r2_barcodes_whitelist", help=_help)
    _help = "Whether to overwrite exitsing files."
    filter_cmd.add_argument("--overwrite", dest="overwrite", action="store_true", help=_help)
    return parser


def main(cli=None):
    print(ct() + "scifi-RNA-seq pipeline")
    args = build_cli().parse_args(cli)
    args.samples = args.samples.split(",") if args.samples != "" else []
    print(args)

    config = get_config(args.config)
    if args.config is not None:
        config = yaml.safe_load(args.config)
    df = pd.read_csv(args.sample_annotation)

    if args.toggle and "toggle" in df.columns:
        df = df.query("toggle == 1")
    if args.samples:
        df = df.loc[df['sample_name'].isin(args.samples), :]

    s = '\n\t - ' + '\n\t - '.join(df['sample_name'])
    print(ct() + f"Samples to submit:{s}")

    for i, sample in df.iterrows():
        print(ct() + f"Doing sample {sample['sample_name']}")

        sample_name = sample['sample_name']
        r1_annotation = pd.read_csv(sample['annotation']).set_index("sample_name")
        sample_out_dir = os.path.join(args.root_output_dir, sample_name)

        if args.command in ["all", "map"]:
            map_command(args, sample_name, sample_out_dir, r1_annotation, config)
        if args.command in ["all", "tracks"]:
            tracks_command(args)
        if args.command in ["all", "filter"]:
            filter_command(
                args, config,
                sample_name, sample_out_dir,
                r1_annotation,
                r1_annotation_file=sample['annotation'],
                r1_attributes=sample['attributes'],
                species_mixture=bool(sample['species_mixing']),
                expected_cell_number=sample['expected_cell_number'],
                correct_r2_barcodes=False, correct_r2_barcodes_file=None)
        if args.command in ["all", "join"]:
            join_command(args)
        if args.command in ["all", "report"]:
            report_command(args)


def get_config(args=None):

    # Config
    # # first load default (distributed with package - to establish structure)
    config = yaml.safe_load()
    # # then load ~/.scifiRNA-seq.config.yaml if present
    if os.path.exists(os.path.expanduser("~/.scifiRNA-seq.config.yaml")):
        c = yaml.safe_load("~/.scifiRNA-seq.config.yaml")
        config.update(c)
    # # then load CLI proved one
    if args:
        if args.config:
            c = yaml.safe_load(args.config)
            config.update(c)
    return config
