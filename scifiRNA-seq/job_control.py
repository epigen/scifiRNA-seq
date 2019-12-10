#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Utilities for job submission and management by the scifi pipeline.
"""


import subprocess


def job_shebang():
    return '#!/bin/env bash\n\ndate\n'


def print_parameters_during_job(job_params):
    return "\n".join([f"#{k} = {v}" for k, v in job_params.items()]) + "\n"


def slurm_echo_array_task_id():
    return 'echo SLURM_ARRAY_TASK_ID = $SLURM_ARRAY_TASK_ID\n'


def slurm_get_array_params_from_array_list(array_file):
    # Get respective line of input based on the $SLURM_ARRAY_TASK_ID var
    txt = (
        f"ARRAY_FILE={array_file}" + "\n"
        "readarray -t ARR < $ARRAY_FILE" + "\n"
        'IFS=" " read -r -a F <<< ${ARR[$SLURM_ARRAY_TASK_ID]}' + "\n"
        "PREFIX=${F[0]}" + "\n"
        "INPUT_BAM=${F[1]}" + "\n")
    return txt


def job_end():
    return '\n\ndate\n'


def write_job_to_file(job, job_file):
    with open(job_file, "w") as handle:
        handle.write(job)


def submit_job(job_file, params, array=None, cmd="sbatch"):
    if array is not None:
        array = f"--array={array} -N 1\\\n"
    params.update({"job_file": job_file, "cmd": cmd})
    cmd = """{cmd} -J {job_name} \\
    -o {log} --time {time} \\
    -c {cpus} --mem {mem} -p {queue} \\
    {array} {job_file}""".format(params)

    subprocess.Popen(cmd.split(" "))


def capture_slurm_job():
    pass
