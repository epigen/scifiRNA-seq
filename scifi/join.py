#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
The main command and supporting functions for the join step of scifi pipeline
"""

import os

from scifi.job_control import (
    job_shebang,
    print_parameters_during_job,
    job_end,
    write_job_to_file,
    submit_job,
)


def join_command(
    sample_name: str,
    sample_out_dir: str,
    r1_attributes: list,
    species_mixture: bool,
    correct_r2_barcodes: bool,
) -> int:
    join_params = dict(cpus=1, mem=8000, queue="shortq", time="00:30:00")

    job_name = f"scifi_pipeline.{sample_name}.join"
    job = os.path.join(sample_out_dir, job_name + ".sh")
    log = os.path.join(sample_out_dir, job_name + ".log")
    params = dict(join_params, job_file=job, log_file=log)

    cmd = job_shebang()
    cmd += print_parameters_during_job(params)
    cmd += join_metrics(
        sample_name=sample_name,
        directory=sample_out_dir,
        r1_attributes=r1_attributes,
        species_mixture=species_mixture,
        exon=False,
        correct_r2_barcodes=correct_r2_barcodes,
    )
    cmd += join_metrics(
        sample_name=sample_name,
        directory=sample_out_dir,
        r1_attributes=r1_attributes,
        species_mixture=species_mixture,
        exon=True,
        correct_r2_barcodes=correct_r2_barcodes,
    )
    cmd += join_expression(
        sample_name=sample_name,
        directory=sample_out_dir,
        r1_attributes=r1_attributes,
        exon=False,
        correct_r2_barcodes=correct_r2_barcodes,
    )
    cmd += join_expression(
        sample_name=sample_name,
        directory=sample_out_dir,
        r1_attributes=r1_attributes,
        exon=True,
        correct_r2_barcodes=correct_r2_barcodes,
    )
    cmd += job_end()
    write_job_to_file(cmd, job)
    submit_job(job, params)
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
        "PREFIX=${F[0]}" + "\n"
        "INPUT_BAM=${F[1]}" + "\n"
    )
    return txt


def join_metrics(
    sample_name: str,
    directory: str,
    r1_attributes: list,
    species_mixture: bool = False,
    exon: bool = False,
    correct_r2_barcodes: bool = False,
) -> str:
    """
    """
    prefix = os.path.join(directory, sample_name)
    filter_exon = "! -name '*exon*'" if not exon else ""
    exon_str = ".exon" if exon else ""
    suffix = "_corrected" if correct_r2_barcodes else ""
    brackets = "{}"

    fields = "r2,read,unique_umi,umi,gene"
    if species_mixture:
        mix_fields = (
            "human,mouse,total,max,ratio,sp_ratio,doublet,unique_fraction,"
            "human_norm,total_norm,max_norm,ratio_norm,sp_ratio_norm"
            "doublet_norm"
        ).split(",")
        mix_fields = (
            "human,mouse,total,max,ratio,sp_ratio,doublet,"
            "read_human,read_mouse,read_total,read_max,read_ratio,read_sp_ratio"
            ",read_doublet,unique_fraction,human_norm,total_norm,max_norm,"
            "ratio_norm,sp_ratio_norm,doublet_norm,read_human_norm,"
            "read_total_norm,read_max_norm,read_ratio_norm,"
            "read_sp_ratio_norm,read_doublet_norm"
        ).split(",")
    else:
        mix_fields = ["unique_fraction"]

    fields = ",".join(fields.split(",") + mix_fields + r1_attributes)

    return f"""
    find {directory}  -mindepth 2 -name '*{exon_str}.metrics{suffix}.csv.gz' {filter_exon} -exec cat {brackets} \\; > {prefix}{exon_str}.metrics{suffix}.csv.gz
    echo '{fields}' > {sample_name}_header1{exon_str}
    gzip {sample_name}_header1{exon_str}
    cat {sample_name}_header1{exon_str}.gz {prefix}{exon_str}.metrics{suffix}.csv.gz > {sample_name}_tmp1{exon_str}
    mv {sample_name}_tmp1{exon_str} {prefix}{exon_str}.metrics{suffix}.csv.gz
    rm {sample_name}_header1{exon_str}.gz\n\n"""


def join_expression(
    sample_name: str,
    directory: str,
    r1_attributes: list,
    exon: bool = False,
    correct_r2_barcodes: bool = False,
) -> str:
    """
    """
    prefix = os.path.join(directory, sample_name)
    filter_exon = "! -name '*exon*'" if not exon else ""
    exon_str = ".exon" if exon else ""
    suffix = "_corrected" if correct_r2_barcodes else ""
    brackets = "{}"

    fields = "r2,gene,umi".split(",") + r1_attributes

    return f"""
    find {directory}  -mindepth 2 -name '*{exon_str}.expression{suffix}.csv.gz' {filter_exon} -exec cat {brackets} \\; > {prefix}{exon_str}.expression{suffix}.csv.gz
    echo '{fields}' > {sample_name}_header2{exon_str}
    gzip {sample_name}_header2{exon_str}
    cat {sample_name}_header2{exon_str}.gz {prefix}{exon_str}.expression{suffix}.csv.gz > {sample_name}_tmp2{exon_str}
    mv {sample_name}_tmp2{exon_str} {prefix}{exon_str}.expression{suffix}.csv.gz
    rm {sample_name}_header2{exon_str}.gz\n\n"""
