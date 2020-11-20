#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
The main command and supporting functions for the filter step of scifi pipeline
"""

import os
from os.path import join as pjoin
import sys
import argparse

import pandas as pd
from scifi import _LOGGER, _CONFIG
from scifi.job_control import (
    job_shebang,
    print_parameters_during_job,
    slurm_echo_array_task_id,
    job_end,
    write_job_to_file,
    submit_job,
)


def filter_command(
    args: argparse.Namespace,
    sample_name: str,
    sample_out_dir: str,
    r1_annotation: pd.DataFrame,
    r1_annotation_file: str,
    r1_attributes: list,
    species_mixture: bool = False,
    expected_cell_number: int = 200000,
    correct_r2_barcodes: bool = False,
    correct_r2_barcodes_file: str = None,
    dry_run: bool = False,
) -> int:
    _LOGGER.debug(f"Running filter command for sample '{sample_name}'")
    filter_params = _CONFIG["resources"]["filter"]
    if correct_r2_barcodes and correct_r2_barcodes_file is None:
        correct_r2_barcodes_file = pjoin(
            sample_out_dir, sample_name + ".fixed_barcodes.mapping.tsv"
        )

    r1_names = list()
    r1_dirs = list()
    for r1_name, r1 in r1_annotation.iterrows():
        r1["sample_name"] = r1.name
        out_dir = pjoin(args.root_output_dir, sample_name, r1_name)
        out_prefix = pjoin(out_dir, r1_name) + ".ALL"
        output_suffix = (
            "metrics" if not args.correct_r2_barcodes else "metrics_corrected"
        )
        out_file = f"{out_prefix}.{output_suffix}.csv.gz"
        _LOGGER.debug(f"Sample '{r1_name}': '{out_prefix}\n{out_file}'")

        if os.path.exists(out_file) and not args.overwrite:
            continue
        r1_names.append(r1_name)
        r1_dirs.append(out_dir)

    if (not r1_names) or (not r1_dirs):
        _LOGGER.debug("Either 'r1_names' or 'r1_dirs' is an empty list.")
        _LOGGER.error(
            "Nothing to process! Output files already found and 'overwrite' option is off!"
        )
        return 1

    if not args.arrayed:
        for r1_name, sample_dir in zip(r1_names, r1_dirs):
            job_name = f"scifi_pipeline.{sample_name}.filter.{r1_name}"
            job = pjoin(sample_out_dir, job_name + ".sh")
            log = pjoin(sample_out_dir, job_name + ".log")
            out_prefix = pjoin(sample_dir, r1_name) + ".ALL"
            params = dict(
                filter_params, job_name=job_name, job_file=job, log_file=log
            )

            cmd = job_shebang()
            cmd += print_parameters_during_job(params)
            cmd += filter_cmd(
                r1_annotation_file,
                r1_attributes,
                prefix=out_prefix,
                sample_name=r1_name,
                exon=False,
                min_umis=_CONFIG["min_umi_output"],
                expected_cell_number=expected_cell_number,
                cell_barcodes="r2",
                species_mixture=species_mixture,
                correct_r2_barcodes=correct_r2_barcodes,
                correct_r2_barcodes_file=correct_r2_barcodes_file,
                overwrite=args.overwrite,
            )
            cmd += filter_cmd(
                r1_annotation_file,
                r1_attributes,
                prefix=out_prefix,
                sample_name=r1_name,
                exon=True,
                min_umis=_CONFIG["min_umi_output"],
                expected_cell_number=expected_cell_number,
                cell_barcodes="r2",
                species_mixture=species_mixture,
                correct_r2_barcodes=correct_r2_barcodes,
                correct_r2_barcodes_file=correct_r2_barcodes_file,
                overwrite=args.overwrite,
            )
            cmd += job_end()
            write_job_to_file(cmd, job)
            submit_job(job, params, dry=args.dry_run)
    else:
        # Write prefix and BAM files to array file
        array_file = pjoin(
            args.root_output_dir,
            sample_name,
            f"scifi_pipeline.{sample_name}.filter.array_file.txt",
        )
        write_array_params(zip(r1_names, r1_dirs), array_file)

        # Now submit job array in chunks of size ``array.size``
        for i in range(0, args.array_size, len(r1_names)):
            array = f"{i}-{i + args.array_size - 1}"
            job_name = f"scifi_pipeline.{sample_name}.filter.{array}"
            job = pjoin(sample_out_dir, job_name + ".sh")
            log = pjoin(sample_out_dir, job_name + ".%a.log")
            params = dict(
                filter_params,
                job_name=job_name,
                job_file=job,
                log_file=log,
                array=array,
            )

            cmd = job_shebang()
            cmd += slurm_echo_array_task_id()
            cmd += get_array_params_from_array_list(array_file)
            cmd += print_parameters_during_job(params)
            cmd += filter_cmd(
                r1_annotation_file,
                r1_attributes,
                prefix=None,
                sample_name=None,
                exon=False,
                min_umis=_CONFIG["min_umi_output"],
                expected_cell_number=expected_cell_number,
                cell_barcodes="r2",
                species_mixture=species_mixture,
                correct_r2_barcodes=correct_r2_barcodes,
                correct_r2_barcodes_file=correct_r2_barcodes_file,
                overwrite=args.overwrite,
            )
            cmd += filter_cmd(
                r1_annotation_file,
                r1_attributes,
                prefix=None,
                sample_name=None,
                exon=True,
                min_umis=_CONFIG["min_umi_output"],
                expected_cell_number=expected_cell_number,
                cell_barcodes="r2",
                species_mixture=species_mixture,
                correct_r2_barcodes=correct_r2_barcodes,
                correct_r2_barcodes_file=correct_r2_barcodes_file,
                overwrite=args.overwrite,
            )
            cmd += job_end()
            write_job_to_file(cmd, job)
            submit_job(job, params, dry=args.dry_run)
    return 0


def write_array_params(params, array_file):
    with open(array_file, "w") as handle:
        for param in params:
            handle.writelines(" ".join(param))


def get_array_params_from_array_list(array_file):
    # Get respective line of input based on the $SLURM_ARRAY_TASK_ID var
    txt = (
        f"ARRAY_FILE={array_file}" + "\n"
        "readarray -t ARR < $ARRAY_FILE" + "\n"
        'IFS=" " read -r -a F <<< ${ARR[$SLURM_ARRAY_TASK_ID]}' + "\n"
        "SAMPLE_NAME=${F[0]}" + "\n"
        "SAMPLE_DIR=${F[1]}" + "\n"
    )
    return txt


def filter_cmd(
    r1_annotation,
    r1_attributes: list,
    prefix=None,
    sample_name=None,
    exon=False,
    min_umis=3,
    expected_cell_number=200000,
    cell_barcodes="r2",
    species_mixture=False,
    correct_r2_barcodes=False,
    correct_r2_barcodes_file=None,
    overwrite=False,
):
    """
    """
    additional_args = ""
    if species_mixture:
        additional_args += "--species-mixture "
    if correct_r2_barcodes:
        additional_args += "--correct-r2-barcodes "
        additional_args += (
            f"--correct-r2-barcode-file {correct_r2_barcodes_file} "
        )
    # align with STAR >=2.7.0e
    if prefix is None:
        prefix = "${SAMPLE_DIR}/${SAMPLE_NAME}"
    if sample_name is None:
        sample_name = "${SAMPLE_NAME}"
    if exon:
        exon = ".exon"
    txt = f"""{sys.executable} -u -m scifi.scripts.summarizer \\
    --r1-annot {r1_annotation} \\
    --r1-attributes {",".join(r1_attributes)} \\
    --cell-barcodes {cell_barcodes} \\
    --only-summary \\
    --no-save-intermediate \\
    --min-umi-output {min_umis} \\
    --expected-cell-number {expected_cell_number} \\
    --no-output-header \\
    --save-gene-expression \\
    --correct-r1-barcodes \\
    {additional_args}--sample-name {sample_name} \\
    {prefix}.*.STAR.Aligned.out{exon}.bam.featureCounts.bam \\
    {prefix}{exon}"""
    return txt + "\n"
