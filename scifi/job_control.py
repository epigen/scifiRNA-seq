#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Utilities for job submission and management by the scifi pipeline.
"""


import subprocess


def job_shebang() -> str:
    return "#!/bin/env bash\n\ndate\n"


def print_parameters_during_job(job_params) -> str:
    return (
        "\n"
        + "\n".join([f"# {k} = {v}" for k, v in job_params.items()])
        + "\n\n"
    )


def slurm_echo_array_task_id() -> str:
    return "echo SLURM_ARRAY_TASK_ID = $SLURM_ARRAY_TASK_ID\n"


def job_end() -> str:
    return "\n\ndate\n\n"


def write_job_to_file(job, job_file) -> None:
    with open(job_file, "w") as handle:
        handle.write(job)


def submit_job(job_file, params, array=None, cmd=None, dry=False) -> None:
    if dry:
        return
    if cmd is None:
        cmd = "sh -e"

    if cmd == "sbatch":
        if array is not None:
            array = f"--array={array} -N 1\\\n"
        params.update({"job_file": job_file, "cmd": cmd})
        cmd = """{cmd} -J {job_name} \\
        -o {log_file} --time {time} \\
        -c {cpus} --mem {mem} -p {queue} \\
        {array}{job_file}""".format(
            array="" if array is not None else "", **params
        )
        subprocess.Popen(cmd.split(" "))
    else:
        subprocess.Popen(cmd.split(" ") + [job_file])


# def capture_slurm_job():
#     pass
