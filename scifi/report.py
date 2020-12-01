#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
The main command and supporting functions for the report step of scifi pipeline
"""

import sys
import os

from scifi.job_control import (
    job_shebang,
    print_parameters_during_job,
    job_end,
    write_job_to_file,
    submit_job,
)


def report_command(
    sample_name: str,
    sample_out_dir: str,
    sample_results_dir: str,
    r1_attributes: list,
    species_mixture: bool,
    correct_r2_barcodes: bool,
):
    report_params = dict(cpus=4, mem=80000, queue="shortq", time="3-00:00:00")

    job_name = f"scifi_pipeline.{sample_name}.report"
    job = os.path.join(sample_out_dir, job_name + ".sh")
    log = os.path.join(sample_out_dir, job_name + ".log")
    params = dict(report_params, job_file=job, log_file=log)

    cmd = job_shebang()
    cmd += print_parameters_during_job(params)
    cmd += report_cmd(
        input_directory=sample_out_dir,
        output_directory=sample_out_dir,
        r1_attributes=r1_attributes,
        species_mixture=species_mixture,
        exon=False,
        correct_r2_barcodes=False,
    )
    cmd += report_cmd(
        input_directory=sample_out_dir,
        output_directory=sample_out_dir,
        r1_attributes=r1_attributes,
        species_mixture=species_mixture,
        exon=True,
        correct_r2_barcodes=False,
    )
    cmd += report_cmd(
        input_directory=sample_out_dir,
        output_directory=sample_out_dir,
        r1_attributes=r1_attributes,
        species_mixture=species_mixture,
        exon=False,
        correct_r2_barcodes=True,
    )
    cmd += report_cmd(
        input_directory=sample_out_dir,
        output_directory=sample_out_dir,
        r1_attributes=r1_attributes,
        species_mixture=species_mixture,
        exon=True,
        correct_r2_barcodes=True,
    )
    cmd += job_end()
    write_job_to_file(cmd, job)
    submit_job(job, params, cmd=args.cmd, dry=args.dry_run)


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
        "PREFIX=${F[0]}" + "\n"
        "INPUT_BAM=${F[1]}" + "\n"
    )
    return txt


def report_cmd(
    sample_name: str,
    input_directory: str,
    output_directory: str,
    r1_attributes: list,
    species_mixture: bool = False,
    exon: bool = False,
    correct_r2_barcodes: bool = False,
) -> str:
    """
    """
    prefix = os.path.join(input_directory, sample_name)
    output_prefix = os.path.join(output_directory, sample_name)
    exon = ".exon" if exon else ""
    suffix = "_corrected" if correct_r2_barcodes else ""
    suffix2 = "corrected." if correct_r2_barcodes else ""
    species_mixture = "--species-mixture" if species_mixture else ""

    return f"""{sys.executable} -u -m scifi.scripts.report \
    {prefix}{exon}.metrics{suffix}.csv.gz \
    {output_prefix}{exon}.{suffix2} \
    --plotting-attributes {",".join(r1_attributes)} {species_mixture}\n\n"""
