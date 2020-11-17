#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
The main command and supporting functions for the mapping step of scifi pipeline
"""

import os
import argparse
from glob import glob
from textwrap import dedent

import pandas as pd

from scifi import _LOGGER
from scifi.job_control import (
    job_shebang,
    print_parameters_during_job,
    slurm_echo_array_task_id,
    job_end,
    write_job_to_file,
    submit_job,
)


def map_command(
    args: argparse.Namespace,
    sample_name: str,
    sample_out_dir: str,
    r1_annotation: pd.DataFrame,
    config: dict,
):
    _LOGGER.debug(f"Running map command for sample '{sample_name}'")
    map_params = dict(cpus=4, mem=60000, queue="shortq", time="08:00:00")

    prefixes = list()
    bams = list()
    _LOGGER.debug("Getting input BAM files for each r1 barcode.")
    attrs = pd.Series(args.input_bam_glob).str.extractall("{(.*?)}").squeeze()
    attrs = set([attrs] if isinstance(attrs, str) else attrs)
    _LOGGER.debug(f"Attributes to use in input BAM files glob: '{attrs}'")
    for r1_name, r1 in r1_annotation.iterrows():
        _LOGGER.debug(f"Getting input BAM files for '{r1_name}'")
        r1["sample_name"] = r1.name
        out_dir = os.path.join(args.root_output_dir, sample_name, r1_name)
        out_prefix = os.path.join(out_dir, r1_name) + ".ALL"
        _LOGGER.debug(f"Prefix for sample '{r1_name}': '{out_prefix}'")

        # get input BAM files
        to_fmt = {attr: r1[attr] for attr in attrs}
        _LOGGER.debug(f"Formatting variables for sample '{r1_name}': '{to_fmt}'")
        bam_file_glob = args.input_bam_glob.format(**to_fmt)
        _LOGGER.debug(f"Glob for BAM files for sample '{r1_name}': '{bam_file_glob}'")
        bam_files = ",".join(glob(bam_file_glob))
        _LOGGER.debug(f"BAM files of sample '{r1_name}': '{bam_files}'")
        prefixes.append(out_prefix)
        bams.append(bam_files)

    if (not prefixes) or (not bam_files):
        _LOGGER.debug("Either 'prefixes' or 'bam_files' is an empty list.")
        _LOGGER.error("Nothing to process! Likely no BAM files were found!")
        return 1

    if not args.arrayed:
        for out_prefix, bam_files in zip(prefixes, bam_files):
            job_name = f"scifi_pipeline.{sample_name}.map.{r1_name}"
            job = os.path.join(sample_out_dir, job_name + ".sh")
            log = os.path.join(sample_out_dir, job_name + ".log")
            params = dict(map_params, job_file=job, log_file=log)

            cmd = job_shebang()
            cmd += print_parameters_during_job(params)
            cmd += star_cmd(
                prefix=out_prefix,
                input_bams=bam_files,
                star_genome_dir=config["star_genome_dir"],
                cpus=4,
                star_exe=config["star_exe"],
            )
            cmd += feature_counts_cmd(gtf_file=config['gtf_file'], prefix=out_prefix, cpus=4, exon=False)
            cmd += link_mapped_file_for_exonic_quantification(prefix=out_prefix)
            cmd += feature_counts_cmd(prefix=out_prefix, gtf_file=config['gtf_file'], cpus=4, exon=True)
            cmd += job_end()
            write_job_to_file(cmd, job)
            submit_job(job, params)
    else:
        # Write prefix and BAM files to array file
        array_file = os.path.join(
            args.root_output_dir, sample_name, f"scifi_pipeline.{sample_name}.map.array_file.txt"
        )
        write_array_params(zip(prefixes, bam_files), array_file)

        # Now submit job array in chunks of size ``array.size``
        for i in range(0, args.array_size, len(bams)):
            array = f"{i}-{i + args.array_size - 1}"
            job_name = f"scifi_pipeline.{sample_name}.map.{array}"
            job = os.path.join(sample_out_dir, job_name + ".sh")
            log = os.path.join(sample_out_dir, job_name + ".%a.log")
            params = dict(map_params, job_file=job, log_file=log, array=array)

            cmd = job_shebang()
            cmd += slurm_echo_array_task_id()
            cmd += get_array_params_from_array_list(array_file)
            cmd += print_parameters_during_job(params)
            cmd += star_cmd(
                prefix=None,
                input_bams=None,
                star_genome_dir=config["star_genome_dir"],
                cpus=4,
                star_exe=config["star_exe"],
            )
            cmd += feature_counts_cmd(prefix=out_prefix, gtf_file=config['gtf_file'], cpus=4, exon=False)
            cmd += link_mapped_file_for_exonic_quantification(prefix=out_prefix)
            cmd += feature_counts_cmd(gtf_file=config['gtf_file'], prefix=None, cpus=4, exon=True)
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


def star_cmd(prefix=None, input_bams=None, star_genome_dir=None, cpus=4, star_exe=None):
    """
    """
    # align with STAR >=2.7.0e
    if prefix is None:
        prefix = "${PREFIX}"
    if input_bams is None:
        input_bams = "${INPUT_BAM}"
    if star_exe is None:
        star_exe = "STAR"
    txt = f"""
    {star_exe} \\
    --runThreadN {cpus} \\
    --genomeDir {star_genome_dir} \\
    --clip3pAdapterSeq AAAAAA \\
    --outSAMprimaryFlag AllBestScore \\
    --outSAMattributes All \\
    --outFilterScoreMinOverLread 0 \\
    --outFilterMatchNminOverLread 0 --outFilterMatchNmin 0 \\
    --outSAMunmapped Within \\
    --outSAMtype BAM Unsorted \\
    --readFilesType SAM SE \\
    --readFilesCommand samtools view -h \\
    --outFileNamePrefix {prefix}.STAR. \\
    --readFilesIn {input_bams}"""
    return dedent(txt) + "\n"


def link_mapped_file_for_exonic_quantification(prefix=None):
    if prefix is None:
        prefix = "${PREFIX}"
    return f"ln -s {prefix}.STAR.Aligned.out.bam \
    {prefix}.STAR.Aligned.out.exon.bam"


def feature_counts_cmd(gtf_file, prefix=None, cpus=4, exon=False):
    if prefix is None:
        prefix = "${PREFIX}"
    # count all reads overlapping a gene
    quant = "exon" if exon else "gene"
    exon = "exon." if exon else ""
    txt = f"""
    featureCounts \\
    -T {cpus} \\
    -F GTF \\
    -t {quant} \\
    -g gene_id \\
    --extraAttributes gene_name \\
    -Q 30 \\
    -s 0 \\
    -R BAM \\
    -a {gtf_file} \\
    -o {prefix}.STAR.featureCounts.quant_gene.{exon}tsv \\
    {prefix}.STAR.Aligned.out.bam"""
    return dedent(txt) + "\n"
