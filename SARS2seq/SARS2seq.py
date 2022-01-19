# pylint: disable=C0103


"""
#placeholder block
"""

import argparse
import multiprocessing
import os
import pathlib
import sys

import snakemake
import yaml

from .functions import MyHelpFormatter, color
from .runconfigs import LoadConf, WriteConfigs
from .samplesheet import WriteSampleSheet
from .update import update
from .userprofile import ReadConfig
from .validatefasta import IsValidFasta
from .version import __version__

yaml.warnings({"YAMLLoadWarning": False})


def get_args(givenargs):
    """
    Parse the commandline args
    """
    # pylint: disable=C0301

    def fasta_input(choices, fname):
        if fname == "NONE":
            return fname
        if os.path.isfile(fname):
            ext = "".join(pathlib.Path(fname).suffixes)
            if ext not in choices:
                raise argparse.ArgumentTypeError(
                    f"Input file doesn't end with one of {choices}"
                )
            return fname
        print(f'"{fname}" is not a file. Exiting...')
        sys.exit(-1)

    def dir_path(arginput):
        if os.path.isdir(arginput):
            return arginput
        print(f'"{arginput}" is not a directory. Exiting...')
        sys.exit(1)

    def currentpath():
        return os.getcwd()

    arg = argparse.ArgumentParser(
        prog="SARS2seq",
        usage="%(prog)s [required arguments] [optional arguments]",
        description="SARS2seq: a dedicated pipeline for analysing SARS-CoV-2 sequencing data in order to generate a consensus sequence specifically tuned to the SARS-CoV-2 virus.",
        formatter_class=MyHelpFormatter,
        add_help=False,
    )

    required_args = arg.add_argument_group("Required arguments")
    optional_args = arg.add_argument_group("Optional arguments")

    required_args.add_argument(
        "--input",
        "-i",
        type=dir_path,
        metavar="DIR",
        help="The input directory with raw fastq(.gz) files",
        required=True,
    )

    required_args.add_argument(
        "--output",
        "-o",
        metavar="DIR",
        type=str,
        default=currentpath(),
        help="Output directory",
        required=True,
    )

    required_args.add_argument(
        "--primers",
        "-pr",
        type=lambda s: fasta_input((".fasta", ".fa"), s),
        metavar="File",
        help="Used primer sequences in FASTA format",
        required=True,
    )

    required_args.add_argument(
        "--platform",
        default="nanopore",
        const="nanopore",
        nargs="?",
        choices=("nanopore", "illumina", "iontorrent"),
        help="Define the sequencing platform that was used to generate the dataset, either being 'nanopore', 'illumina' or 'iontorrent', see the docs for more info",
        required=True,
    )

    required_args.add_argument(
        "--amplicon-type",
        "-at",
        default="end-to-end",
        const="end-to-end",
        nargs="?",
        choices=("end-to-end", "end-to-mid"),
        help="Define the amplicon-type, either being 'end-to-end' or 'end-to-mid', see the docs for more info",
        required=True,
    )

    optional_args.add_argument(
        "--threads",
        "-t",
        default=min(multiprocessing.cpu_count(), 128),
        metavar="N",
        type=int,
        help=f"Number of local threads that are available to use.\nDefault is the number of available threads in your system ({min(multiprocessing.cpu_count(), 128)})",
    )

    optional_args.add_argument(
        "--version",
        "-v",
        version=__version__,
        action="version",
        help="Show the SARS2seq version and exit",
    )

    optional_args.add_argument(
        "--help",
        "-h",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit",
    )

    optional_args.add_argument(
        "--dryrun",
        action="store_true",
        help="Run the workflow without actually doing anything",
    )

    optional_args.add_argument(
        "--skip-updates", action="store_true", help="Skip the update check",
    )

    if len(givenargs) < 1:
        print(
            f"{arg.prog} was called but no arguments were given, please try again\n\tUse '{arg.prog} -h' to see the help document"
        )
        sys.exit(1)
    else:
        flags = arg.parse_args(givenargs)

    return flags


def CheckInputFiles(indir):
    """
    Check if the input files are valid fastq files
    """
    allowedextensions = [".fastq", ".fq", ".fastq.gz", ".fq.gz"]
    foundfiles = []

    for filenames in os.listdir(indir):
        extensions = "".join(pathlib.Path(filenames).suffixes)
        foundfiles.append(extensions)

    return bool(any(i in allowedextensions for i in foundfiles))


def main():
    """
    SARS2seq starting point
    --> Fetch and parse arguments
    --> check validity
    --> Read (or write, if necessary) the user-config files
    --> Change working directories and make necessary local files for snakemake
    --> Run snakemake with appropriate settings
    """

    ##> Check the default userprofile, make it if it doesn't exist
    conf = ReadConfig(os.path.expanduser("~/.SARS2seq_defaultprofile.ini"))

    flags = get_args(sys.argv[1:])

    if not flags.skip_updates:
        update(sys.argv, conf)

    inpath = os.path.abspath(flags.input)
    # refpath = os.path.abspath(flags.reference)
    if flags.primers != "NONE":
        primpath = os.path.abspath(flags.primers)
    else:
        primpath = "NONE"
    outpath = os.path.abspath(flags.output)

    here = os.path.abspath(os.path.dirname(__file__))

    Snakefile = os.path.join(here, "workflow", "workflow.smk")

    ##@ check if the input directory contains valid files
    if CheckInputFiles(inpath) is False:
        print(
            f"""
{color.RED + color.BOLD}"{inpath}" does not contain any valid FastQ files.{color.END}
Please check the input directory. Exiting...
            """
        )
        sys.exit(-1)
    else:
        print(
            f"""
{color.GREEN}Valid input files were found in the input directory{color.END} ('{inpath}')
            """
        )

    if IsValidFasta(primpath) is False:
        print(
            f"""
{color.RED + color.BOLD}
The given fasta with primer sequences contains illegal characters in its sequences.
{color.END}
Please check the primer fasta and try again. Exiting...
            """
        )
        sys.exit(1)

    ##@ check if the output dir exists, create if not
    ##@ change the working directory
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    # copy_tree(os.path.join(here, 'envs'), os.path.join(outpath, 'envs'))

    if not os.getcwd() == outpath:
        os.chdir(outpath)
    workdir = outpath

    samplesheet = WriteSampleSheet(inpath, flags.platform)
    snakeparams, snakeconfig = WriteConfigs(
        conf,
        flags.threads,
        os.getcwd(),
        flags.platform,
        primpath,
        samplesheet,
        flags.amplicon_type,
        flags.dryrun,
    )

    parsedconfig = LoadConf(snakeconfig)

    if conf["COMPUTING"]["compmode"] == "local":
        snakemake.snakemake(
            Snakefile,
            workdir=workdir,
            cores=parsedconfig["cores"],
            use_conda=parsedconfig["use-conda"],
            conda_frontend="mamba",
            jobname=parsedconfig["jobname"],
            latency_wait=parsedconfig["latency-wait"],
            dryrun=parsedconfig["dryrun"],
            configfiles=[snakeparams],
            restart_times=3,
        )
    if conf["COMPUTING"]["compmode"] == "grid":
        snakemake.snakemake(
            Snakefile,
            workdir=workdir,
            cores=parsedconfig["cores"],
            nodes=parsedconfig["cores"],
            use_conda=parsedconfig["use-conda"],
            conda_frontend="mamba",
            jobname=parsedconfig["jobname"],
            latency_wait=parsedconfig["latency-wait"],
            drmaa=parsedconfig["drmaa"],
            drmaa_log_dir=parsedconfig["drmaa-log-dir"],
            dryrun=parsedconfig["dryrun"],
            configfiles=[snakeparams],
            restart_times=3,
        )

    if parsedconfig["dryrun"] is False:
        snakemake.snakemake(
            Snakefile,
            workdir=workdir,
            report="results/snakemake_report.html",
            configfiles=[snakeparams],
            quiet=True,
        )
