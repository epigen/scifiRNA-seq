#!/usr/bin/env python

"""
sciRNA-seq barcode extraction and correction script.
"""

import sys
from argparse import ArgumentParser
import numpy as np
import pandas as pd
import pysam
import time

__author__ = "Andre Rendeiro"
__copyright__ = "Copyright 2018, Andre Rendeiro"
__credits__ = []
__license__ = "GPL2"
__version__ = "0.2"
__maintainer__ = "Andre Rendeiro"
__email__ = "arendeiro@cemm.oeaw.ac.at"
__status__ = "Development"


def main():
    # Parse command-line arguments
    parser = ArgumentParser(
        prog="python scirnaseq.barcodes.py",
        description="\n".join([
            "sciRNA-seq script from Bock lab. " +
            "See https://github.com/epigen/sciRNA-seq for specific documentation."])
    )
    parser = arg_parser(parser)
    # args = parser.parse_args(
    #    ("-a metadata/sciRNA-seq.oligos_2018-09-17.csv " +
    #    "--mode slim "
    #    "/scratch/lab_bsf/samples/BSF_0513_XXXXXXXXX/BSF_0513_XXXXXXXXX_4_samples/BSF_0513_XXXXXXXXX_4#SCI_011_3_gate_more_S43637.bam").split(" "))
    args = parser.parse_args()

    args.barcodes = args.barcodes.split(",")
    args.barcode_tags = args.barcode_tags.split(",")
    args.barcode_lengths = args.barcode_lengths.split(",")
    args.barcode_lengths = [int(x) for x in args.barcode_lengths]
    args.correct_barcodes = args.correct_barcodes.split(",")
    msg = "Barcode tags ({}) and lengths ({}) are not of the same number as the number of barcodes to extract ({})!".format(
        ", ".join(args.barcode_tags), ", ".join([str(x) for x in args.barcode_lengths]), ", ".join(args.barcodes)
    )
    assert len(args.barcodes) == len(args.barcode_tags) == len(args.barcode_lengths), msg
    msg = "Barcodes to correct ({}) are not among barcodes to extract ({})!".format(
        ", ".join(args.correct_barcodes), ", ".join(args.barcodes)
    )
    assert all([x in args.barcodes for x in args.correct_barcodes]), msg

    print("# " + time.asctime() + " - Start.")
    print(args)

    annotation = pd.read_csv(args.annotation)

    # Get barcodes and annotate with mismatches
    cells = annotate_barcodes(
        extract_barcodes(
            args.input_file,
            barcodes=args.barcodes,
            barcode_tags=args.barcode_tags,
            barcode_lengths=args.barcode_lengths,
            start=args.start, end=args.end),
        annotation,
        barcodes=args.correct_barcodes)

    # Slim down output if required
    if args.mode == "slim":
        # remove Ns
        cells = cells.loc[
            ~(cells.loc[:, cells.columns[cells.columns.str.contains('_contains_N')]] == 1)
            .any(axis=1), :]

        if args.max_mismatches > 0:
            for barcode in args.correct_barcodes:
                f = (
                    (cells[barcode + "_mismatches"] > 0) &
                    (cells[barcode + "_closest"] != "X") &
                    (cells[barcode + "_mismatches"] <= args.max_mismatches))
                cells.loc[f, barcode] = cells.loc[f, barcode + "_closest"]
                cells = cells.drop(barcode + "_closest", axis=1)

        # remove not matching any barcode with required mismatch threshold
        cells = cells.loc[
            ~(cells.loc[
                :,
                cells.columns[cells.columns.str.contains('_mismatches')]] > args.max_mismatches)
            .any(axis=1), :]
        # remove unused column
        cells = cells.drop(cells.columns[cells.columns.str.contains('_')], axis=1)

    # Save
    o = ['read'] + args.barcodes
    o += [x for x in cells.columns if x not in o]
    print("# Saving to file: " + args.output_file)
    cells[o].sort_values("read").to_csv(args.output_file, index=False, compression="gzip")
    print("# " + time.asctime() + " - Done.")


def arg_parser(parser):
    """
    Global options.
    """
    parser.add_argument(
        dest="input_file",
        help="Input BAM file with reads to process.",
        type=str)
    parser.add_argument(
        "-a", "--annotation",
        dest="annotation",
        help="CSV file with barcode annotation.",
        type=str)
    default = "sciRNA-seq.barcodes.csv.gz"
    parser.add_argument(
        "-o", "--output",
        dest="output_file",
        help="Output file with barcodes. Default is '{}'.".format(default),
        default=default,
        type=str)
    default = ['round1', 'round2', 'umi']
    parser.add_argument(
        "-b", "--barcodes",
        dest="barcodes",
        help="Existing barcodes in the experiment. " +
             "Comma-separated list of barcodes. " +
             "Default is '{}'".format(",".join(default)),
        default=",".join(default),
        type=str)
    default = ['r1', 'r2', 'RX']
    parser.add_argument(
        "--barcode-tags",
        dest="barcode_tags",
        help="BAM file barcodes tags in the experiment. " +
             "Comma-separated list of barcodes, should match number of values in '--barcodes'. " +
             "Default is '{}'".format(",".join(default)),
        default=",".join(default),
        type=str)
    default = ['11', '16', '8']
    parser.add_argument(
        "--barcode-lengths",
        dest="barcode_lengths",
        help="Length of barcodes tags in the experiment. " +
             "Comma-separated list of barcodes, should match number of values in '--barcodes'. " +
             "Default is '{}'".format(",".join(default)),
        default=",".join(default),
        type=str)
    default = 0
    parser.add_argument(
        "--start",
        dest="start",
        help="Start line of BAM file to begin processing. " +
             "Default is first line (i.e. {}th line).".format(default),
        default=default,
        type=int)
    default = 1e100
    parser.add_argument(
        "--end",
        dest="end",
        help="End line of BAM file to finish processing. " +
             "Default is whole file (i.e. {}th  line)".format(default),
        default=default,
        type=int)
    parser.add_argument(
        "-d", "--dry-run",
        dest="dry_run",
        help="Whether not to actually do any work but just check files.",
        action="store_true")
    choices = ["fat", "slim"]
    parser.add_argument(
        "--mode",
        dest="mode",
        default=choices[1],
        choices=choices,
        help="Whether barcode correction should be applied and no record of original " +
             "barcode should be kept ('slim' mode), or all records should be kept in " +
             "the output file ('fat' mode). " +
             "One of ['{}']. Default is '{}'.".format(
                "', '".join(choices), choices[1]),
        type=str)
    default = 1
    parser.add_argument(
        "--max-mismatches",
        dest="max_mismatches",
        help="Maximum mismatches to allow correction. Default {}".format(default),
        default=default,
        type=int)
    default = ['round1']
    parser.add_argument(
        "--correct-barcodes",
        dest="correct_barcodes",
        help="Cell-labeling barcodes to match against reference and correct. " +
             "Comma-separated list of barcodes matching '--barcodes'. " +
             "Default is '{}'".format(",".join(default)),
        default=",".join(default),
        type=str)

    return parser


def extract_barcodes(
        input_file,
        barcodes=['round1', 'round2', 'umi'],
        barcode_tags=['r1', 'r2', 'RX'],
        barcode_lengths=[11, 16, 8],
        start=0, end=1e100):
    """
    """
    from collections import Counter
    input_handle = pysam.AlignmentFile(input_file, mode="rb", check_sq=False)

    errors = Counter()
    cells = list()
    i = 0
    print("# " + time.asctime() + " - Starting to extract barcodes '{}'.".format(", ".join(barcodes)))
    for read in input_handle:
        i += 1

        if (start > (i - 1)):
            continue
        if (end < (i - 1)):
            break

        if (i - 1) % 1000000 == 0:
            print(i)

        res = [read.qname]
        for barcode, tag, length in zip(barcodes, barcode_tags, barcode_lengths):
            if read.has_tag(tag):
                seq = read.get_tag(tag)
                if len(seq) != length:
                    errors[barcode + "_not_right_length"] += 1
                else:
                    res.append(seq)
            else:
                errors[barcode + "_no_tag"] += 1
        cells.append(res)

    print("# " + time.asctime() + " - Done extracting barcodes.")
    print("## Errors:")
    print(errors)
    return pd.DataFrame(data=cells, columns=["read"] + barcodes)


def annotate_barcodes(cells, annotation, barcodes=["round1", "round2"]):
    def count_mismatches(a, b):
        return sum(a != b for a, b in zip(a, b))

    # fraction mapping to annotation
    print("# " + time.asctime() + " - Starting to annotate barcodes.")
    print("## Matching to annotation.")
    for a, barcode in enumerate(barcodes):
        print(" - " + barcode)
        cells.loc[:, barcode + "_correct"] = (
            cells[barcode].isin(
                annotation.loc[annotation["barcode_type"] == barcode, "barcode_sequence"])).astype(int)
        cells.loc[:, barcode + '_contains_N'] = cells[barcode].str.contains("N").astype(int)

    print("# " + time.asctime() + " - Starting to correct barcodes.")
    for barcode in barcodes[::-1]:
        print("## " + time.asctime() + " - " + barcode)
        ref = annotation.loc[annotation["barcode_type"] == barcode, "barcode_sequence"]
        print(" - Creating query")
        query = cells.loc[
            (cells[barcode + "_correct"] == 0) & (cells[barcode + "_contains_N"] == 0),
            barcode].drop_duplicates()
        print(" - Finding mismatches against reference")
        mis = pd.DataFrame([
            (q, count_mismatches(q, b), b)
            for b in ref
            for q in query], columns=[barcode, barcode + "_mismatches", barcode + '_closest'])
        print(" - " + time.asctime() + " -  Finding minimal match.")
        fix = mis.loc[mis.groupby(barcode)[barcode + "_mismatches"].idxmin()]
        print(" - " + time.asctime() + " -  Merging correct to reference.")
        cells = cells.set_index(barcode).join(fix.set_index(barcode)).reset_index()

        # have all entries with same type (no NAs)
        cells[barcode + "_mismatches"] = cells[barcode + "_mismatches"].fillna(0)
        cells[barcode + "_mismatches"] = cells[barcode + "_mismatches"].astype(int)
        cells[barcode + "_closest"] = cells[barcode + "_closest"].fillna("X")

    print("# Time: " + time.asctime() + " - Finished correcting barcodes.")
    return cells.sort_values("read")


if __name__ == '__main__':
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print("Program canceled by user!")
        sys.exit(1)
