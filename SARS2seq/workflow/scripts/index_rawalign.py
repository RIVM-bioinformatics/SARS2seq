import argparse
import os
import pathlib
import sys

import pandas as pd
import pysam
import pysamstats


def GetArgs(sysargs):
    def checkbam(fname):
        if os.path.isfile(fname):
            ext = "".join(pathlib.Path(fname).suffix)
            if ext != ".bam":
                args.error(f"Input file ({fname}) doesn't seem to be a BAM-file.")
            return fname
        else:
            print(f'"{fname}" is not a file. Exiting...')
            sys.exit(-1)

    args = argparse.ArgumentParser()

    args.add_argument(
        "-i",
        "--input",
        type=lambda s: checkbam(s),
        required=True,
        help="Input bam file",
    )
    args.add_argument(
        "-r", "--reference", type=str, required=True, help="Reference fasta file"
    )
    args.add_argument("-o", "--output", type=str, required=True, help="Output file")

    flags = args.parse_args(sysargs)

    return flags


def Readbam(f):
    return pysam.AlignmentFile(f, "rb")


def BuildIndex(bamfile, ref):
    columns = ["coverage", "A", "T", "C", "G", "X", "I"]
    p_index = pd.DataFrame(columns=columns)

    for r in pysamstats.stat_pileup(
        type="variation",
        alignmentfile=Readbam(bamfile),
        stepper="nofilter",
        fafile=ref,
        pad=True,
        one_based=True,
        max_depth=1000000000,
    ):
        p_index.loc[r["pos"]] = (
            [r["reads_all"]]
            + [r["A"]]
            + [r["T"]]
            + [r["C"]]
            + [r["G"]]
            + [r["deletions"]]
            + [r["insertions"]]
        )

    return p_index


if __name__ == "__main__":
    flags = GetArgs(sys.argv[1:])

    bam_index = BuildIndex(flags.input, flags.reference)

    bam_index.to_csv(flags.output, sep=",", index=True, compression="gzip")
