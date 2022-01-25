import argparse
import os
import pathlib
import sys

import pandas as pd
import pysam


def GetArgs(sysargs):
    def checkbam(fname):
        if os.path.isfile(fname):
            ext = "".join(pathlib.Path(fname).suffix)
            if ext != ".bam":
                args.error(f"Input file ({fname}) doesn't seem to be a BAM-file.")
            return fname
        print(f'"{fname}" is not a file. Exiting...')
        sys.exit(-1)

    args = argparse.ArgumentParser()

    args.add_argument(
        "-i",
        "--input",
        type=checkbam,
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
    bamfile = pysam.AlignmentFile(bamfile, "rb")
    ref_fasta = pysam.FastaFile(ref)
    ref_length = ref_fasta.lengths[0]

    pileup = bamfile.pileup(stepper="nofilter", max_depth=10000000, min_base_quality=0)

    def parse_query_sequences(l):
        coverage = a = c = t = g = x = i = 0
        for b in l:
            coverage += 1
            if b == "*":
                x += 1
            elif b[0].lower() == "a":
                a += 1
            elif b[0].lower() == "t":
                t += 1
            elif b[0].lower() == "c":
                c += 1
            elif b[0].lower() == "g":
                g += 1

            # It is important to count the insertions seperately
            if "+" in b:
                i += 1
        return coverage, a, t, c, g, x, i

    columns = ["pos", "coverage", "A", "T", "C", "G", "X", "I"]

    # 1 Is added to the position because our index starts at 1
    p_index = pd.DataFrame(
        (
            (p.pos + 1,) + parse_query_sequences(p.get_query_sequences(add_indels=True))
            for p in pileup
        ),
        columns=columns,
    )

    # Since the pileup does not return positions without any reads mapped, we have to
    # fill these with zeroes
    missing_positions = set(range(1, ref_length + 1)) - set(p_index.pos)
    missing_p_index = pd.DataFrame(
        ((i, 0, 0, 0, 0, 0, 0, 0) for i in missing_positions), columns=columns
    )
    p_index = p_index.append(missing_p_index).set_index("pos").sort_index()
    p_index.index.name = None

    return p_index


if __name__ == "__main__":
    flags = GetArgs(sys.argv[1:])

    bam_index = BuildIndex(flags.input, flags.reference)

    bam_index.to_csv(flags.output, sep=",", index=True, compression="gzip")
