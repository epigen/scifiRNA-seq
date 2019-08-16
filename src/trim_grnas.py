#! /bin/env python

import sys
from argparse import ArgumentParser
import time

import pysam


def parse_args():
    parser = ArgumentParser()
    parser.add_argument(
        dest="input_bam",
        help="Input BAM file",
    )
    parser.add_argument(
        dest="output_bam",
        help="Output BAM file",
    )
    parser.add_argument(
        "-s", "--sequence",
        dest="sequence",
        help="Output BAM file",
        default="GTGGAAAGGACGAAACACCG"
    )
    # # Example:
    # args = parser.parse_args([
    #     "--sequence", "GTGGAAAGGACGAAACACCG",
    #     "/scratch/users/dbarreca/private/custom_demux/scRNA/BSF_0662_HCWH5DRXX/BSF_0662_HCWH5DRXX_1_samples/BSF_0662_HCWH5DRXX_1#PD188_SingleCell_Ampligase_gRNA.bam",
    #     "/home/arendeiro/projects/sci-rna/data/PD187/PD188_SingleCell_Ampligase_gRNA/PD188_SingleCell_Ampligase_gRNA.input.bam"])
    args = parser.parse_args()

    return args


def main():
    print(f"# {time.asctime()} - Starting!")
    args = parse_args()
    trim_seq(args.input_bam, args.output_bam, args.sequence)
    print(f"# {time.asctime()} - Finished!")


def trim_seq(input_bam, output_bam, sequence):

    input_ = pysam.AlignmentFile(input_bam, check_sq=False)
    output_ = pysam.AlignmentFile(output_bam, "wb", template=input_)
    n = len(sequence)

    for r in input_:
        try:
            i = r.seq.index(sequence)
        except ValueError:
            continue
        q = r.query_qualities[i + n:]
        r.seq = r.seq[i + n:]
        r.query_qualities = q
        output_.write(r)

    input_.close()
    output_.close()


if __name__ == "__main__":
    sys.exit(main())
