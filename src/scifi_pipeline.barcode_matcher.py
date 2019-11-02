#!/usr/bin/env python

"""
A helper script to match barcodes to a reference.
"""


import sys
import os
from argparse import ArgumentParser
import time

import numpy as np
import pandas as pd

import pyximport
if os.path.exists("/home/arendeiro/sci-rna/src"):
    sys.path.insert(0, "/home/arendeiro/sci-rna/src")

pyximport.install()
from scifi_utils_c import edit_distance


def parse_arguments():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--metrics", required=True)
    default = "/home/arendeiro/sci-rna/metadata/737K-cratac-v1.txt"
    parser.add_argument("--output-prefix")
    parser.add_argument("--r2-barcodes", default=default)
    parser.add_argument("--total-chunks", type=int, default=250)
    parser.add_argument("--start-chunk", type=int, required=True)
    parser.add_argument("--batch-size", type=int, required=True)
    parser.add_argument("--max-dist", type=int, default=1)
    parser.add_argument("--parallel", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    # parser.add_argument("--job-array", action="store_true")
    return parser


def closest_match(b, whitelist, max_dist=1):
    m = [edit_distance(b, x) for x in whitelist]
    r = min(m)
    if r <= max_dist:
        return np.argmin(m)
    else:
        return None


def main(cli=None):
    # cli = [
    #     '--prefix', '/home/arendeiro/sci-rna/data/PD193_fivelines_383k/PD193_fivelines_383k.',
    #     '--metrics', '/home/arendeiro/sci-rna/data/PD193_fivelines_383k/PD193_fivelines_383k.metrics.csv.gz',
    #     '--chunks', '10',
    #     '--current-chunk', '0'
    # ]
    # args = parse_arguments().parse_args(cli)
    args = parse_arguments().parse_args()
    if args.output_prefix is None:
        args.output_prefix = args.metrics.replace("metrics.csv.gz", "")
    # if args.job_array:
    #     args.current_chunk = os.env['SLURM_ARRAY_TASK_ID']
    # else:
    #     if args.current_chunk is None:
    #         print("Option --current_chunk is mandatory!")
    #         sys.exit(1)

    if args.parallel:
        import multiprocessing
        import parmap

    print(f"# {time.asctime()} - Reading in input metrics file.")
    metrics = pd.read_csv(args.metrics)
    print(f"# {time.asctime()} - Reading in barcode whitelist.")
    r2_barcodes = pd.read_csv(
        args.r2_barcodes, header=None, squeeze=True).tolist()

    # get non-matching
    match = metrics['r2'].isin(r2_barcodes)
    to_match = metrics.loc[~match].sort_values("umi")
    to_match = to_match.groupby('r2')['umi'].sum().sort_values()

    # convert to bytes
    barcodes = [bytes(x, 'ascii') for x in to_match.index]
    whitelist = [bytes(x, 'ascii') for x in r2_barcodes]

    # get respective chunks
    print(f"# {time.asctime()} - Chunking the input.")
    barcode_chunks = np.array_split(barcodes, args.total_chunks)

    for chunk in range(args.start_chunk, args.start_chunk + args.batch_size):
        c = f"{chunk}/{args.total_chunks}"
        output_file = os.path.join(
            args.output_prefix + f"fixed_barcodes.{args.total_chunks}-{chunk}.tsv")

        if not args.overwrite:
            if os.path.exists(args.output_prefix):
                print(f"# {time.asctime()} - Output file for chunk '{c}' already exits!")
                continue

        print(f"# {time.asctime()} - Doing chunk '{c}'!")

        chunk_barcodes = barcode_chunks[chunk].tolist()
        print(
            f"# {time.asctime()} - Getting matches for chunk"
            f" '{c}' with {len(chunk_barcodes)} barcodes.")

        # match barcodes
        if args.parallel:
            matches = parmap.map(
                closest_match, chunk_barcodes,
                whitelist=whitelist, max_dist=args.max_dist)

            # write out
            print(f"# {time.asctime()} - Writing to file: '{output_file}'")
            with open(output_file, 'wb') as handle:
                for k, v, in zip(chunk_barcodes, matches):
                    v = b"" if not v else whitelist[v]
                    handle.write(k + b"\t" + v + b"\n")
        else:
            print(f"# {time.asctime()} - Running in serial and writing as it goes to: '{output_file}'")
            with open(output_file, 'wb') as handle:
                for b in chunk_barcodes:
                    v = closest_match(b, whitelist=whitelist, max_dist=args.max_dist)
                    v = b"" if not v else whitelist[v]
                    handle.write(b + b"\t" + v + b"\n")

        print(f"# {time.asctime()} - Finished chunk '{c}'!")
    print(f"# {time.asctime()} - Finished all chunks!")


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print(f"# {time.asctime()} - KeyboardInterrupt by user.")
        sys.exit(1)
