"""
#placeholder block
"""

# * base imports
import os
import sys
import argparse
import pathlib
import multiprocessing
import snakemake
import yaml

yaml.warnings({'YAMLLoadWarning': False})

# * relative imports
from .version import __version__
from .functions import MyHelpFormatter, color
from .validatefasta import IsValidFasta
from .samplesheet import WriteSampleSheet
from .userprofile import ReadConfig
from .runconfigs import WriteConfigs


def get_args(givenargs):
    """
    Parse the commandline args
    """

    def fasta_input(choices, fname):
        if os.path.isfile(fname):
            ext = "".join(pathlib.Path(fname).suffixes)
            if ext not in choices:
                arg.error("Input file doesn't end with one of {}".format(choices))
            return fname
        else:
            print(f'"{fname}" is not a file. Exiting...')
            sys.exit(-1)
            
    def dir_path(arginput):
        if os.path.isdir(arginput):
            return arginput
        else:
            print(f'\"{arginput}\" is not a directory. Exiting...')
            sys.exit(1)
    
    def currentpath():
        return os.getcwd()

    arg = argparse.ArgumentParser(
        prog="SARS2seq",
        usage="%(prog)s [required options] [optional arguments]",
        description="SARS2seq: a dedicated pipeline for analysing SARS-CoV-2 sequencing data in order to generate a consensus sequence specifically tuned to the SARS-CoV-2 virus.",
        formatter_class=MyHelpFormatter,
        add_help=False,
    )

    arg.add_argument(
        "--input",
        "-i",
        type=dir_path,
        metavar="DIR",
        help="The input directory with raw fastq(.gz) files",
        required=True,
    )

    arg.add_argument(
        "--output",
        "-o",
        metavar="DIR",
        type=str,
        default=currentpath(),
        help="Output directory",
        required=True,
    )

    arg.add_argument(
        "--reference",
        "-ref",
        type=lambda s: fasta_input((".fasta", ".fa"), s),
        metavar="File",
        help="Input Reference genome in FASTA format",
        required=True,
    )

    arg.add_argument(
        "--primers",
        "-pr",
        type=lambda s: fasta_input((".fasta", ".fa"), s),
        metavar="File",
        help="Used primer sequences in FASTA format",
        required=True,
    )
    
    arg.add_argument(
        "--platform",
        default="nanopore",
        const="nanopore",
        nargs='?',
        choices=('nanopore', 'illumina'),
        help="Define the sequencing platform that was used to generate the dataset, either being 'nanopore' or 'illumina', see the docs for more info",
        required=True
    )

    arg.add_argument(
        "--amplicon-type",
        "-at",
        default="end-to-end",
        const="end-to-end",
        nargs="?",
        choices=("end-to-end", "end-to-mid"),
        help="Define the amplicon-type, either being 'end-to-end' or 'end-to-mid', see the docs for more info",
        required=True,
    )
    
    arg.add_argument(
        "--threads",
        "-t",
        default=min(multiprocessing.cpu_count(), 128),
        metavar="N",
        type=int,
        help=f"Number of local threads that are available to use.\nDefault is the number of available threads in your system ({min(multiprocessing.cpu_count(), 128)})"
    )

    arg.add_argument(
        "--version",
        "-v",
        version=__version__,
        action="version",
        help="Show the SARS2seq version and exit",
    )

    arg.add_argument(
        "--help",
        "-h",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit",
    )
    
    arg.add_argument(
        "--dryrun",
        action='store_true',
        help="Run the workflow without actually doing anything"
    )

    if len(givenargs)<1:
        print(f"{arg.prog} was called but no arguments were given, please try again\n\tUse '{arg.prog} -h' to see the help document")
        sys.exit(1)
    else:
        flags = arg.parse_args(givenargs)

    return flags

def CheckInputFiles(dir):
    allowedextensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz']
    foundfiles = []
    
    for filenames in os.listdir(dir):
        extensions = ''.join(pathlib.Path(filenames).suffixes)
        foundfiles.append(extensions)
        
    if any(i in allowedextensions for i in foundfiles) is True:
        return True
    else:
        return False

def main():
    flags = get_args(sys.argv[1:])
    
    inpath = os.path.abspath(flags.input)
    refpath = os.path.abspath(flags.reference)
    primpath = os.path.abspath(flags.primers)
    outpath = os.path.abspath(flags.output)
    
    here = os.path.abspath(os.path.dirname(__file__))
    
    Snakefile = os.path.join(here, "workflow", "workflow.smk")
    
    ##@ check if the input directory contains valid files
    if CheckInputFiles(inpath) is False:
        print(f'{color.RED + color.BOLD}\"{inpath}\" does not contain any valid FastQ files.{color.END}\nPlease check the input directory. Exiting...')
        sys.exit(-1)
    else:
        print(f'{color.GREEN}Valid input files were found in the input directory{color.END} ({inpath})')
    
    ##> Check the default userprofile, make it if it doesn't exist
    conf = ReadConfig(os.path.expanduser('~/.SARS2seq_defaultprofile.ini'))
    
    ##@ validate the reference and primer files
    if IsValidFasta(refpath) is False:
        print(f'{color.RED + color.BOLD}The given reference fasta contains illegal characters in its sequence.{color.END}\nPlease check the reference fasta and try again. Exiting...')
        sys.exit(1)
    
    if IsValidFasta(primpath) is False:
        print(f'{color.RED + color.BOLD}The given fasta with primer sequences contains illegal characters in its sequences.{color.END}\nPlease check the primer fasta and try again. Exiting...')
        sys.exit(1)
    
    ##@ check if the output dir exists, create if not
    ##@ change the working directory
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    
    #copy_tree(os.path.join(here, 'envs'), os.path.join(outpath, 'envs'))
    
    if not os.getcwd() == outpath:
        os.chdir(outpath)
    workdir = outpath
        
    samplesheet = WriteSampleSheet(inpath, flags.platform)
    snakeparams, snakeconfig = WriteConfigs(conf, flags.threads, os.getcwd(), flags.platform, refpath, primpath, samplesheet, flags.amplicon_type, flags.dryrun)
    
    openedconfig = open(snakeconfig)
    parsedconfig = yaml.load(openedconfig, Loader=yaml.FullLoader)
    
    if conf['COMPUTING']['compmode'] == 'local':
        snakemake.snakemake(Snakefile, 
                            workdir=workdir, 
                            cores=parsedconfig['cores'], 
                            use_conda=parsedconfig['use-conda'],
                            conda_frontend="mamba",
                            jobname=parsedconfig['jobname'],
                            latency_wait=parsedconfig['latency-wait'],
                            dryrun=parsedconfig['dryrun'],
                            configfiles=[snakeparams]
                            )
    if conf['COMPUTING']['compmode'] == 'grid':
        snakemake.snakemake(Snakefile,
                            workdir=workdir,
                            cores=parsedconfig['cores'], 
                            use_conda=parsedconfig['use-conda'],
                            conda_frontend="mamba",
                            jobname=parsedconfig['jobname'],
                            latency_wait=parsedconfig['latency-wait'],
                            drmaa=parsedconfig['drmaa'],
                            drmaa_log_dir=parsedconfig['drmaa-log-dir'],
                            dryrun=parsedconfig['dryrun'],
                            configfiles=[snakeparams])