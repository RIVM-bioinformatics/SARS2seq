# SARS2seq

[![CodeFactor](https://www.codefactor.io/repository/github/rivm-bioinformatics/sars2seq/badge)](https://www.codefactor.io/repository/github/rivm-bioinformatics/sars2seq)
![GitHub top language](https://img.shields.io/github/languages/top/RIVM-bioinformatics/SARS2seq)
![Snakemake](https://img.shields.io/badge/snakemake-6.4.1-brightgreen.svg?style=flat-square)

![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/RIVM-bioinformatics/SARS2seq?include_prereleases)
![GitHub](https://img.shields.io/github/license/RIVM-bioinformatics/SARS2seq)

---

SARS2seq is a pipeline designed to process raw FastQ data from targeted SARS2-CoV-2 sequencing and generate biologically correct consensus sequences of the SARS-CoV-2 genome.

SARS2seq performs high speed data quality control and data cleanup, and high accuracy removal of primer-sequences from NGS reads. As well as alignment of reads and generation of a consensus sequence using a custom consensus caller which accounts for sequencing errors and alignment artefacts.


SARS2seq is able to run both on a standalone (linux) computer, as well as High-Performance Computing (HPC) infrastructures.

---

## Download & installation

***Download and installation instructions through conda will be made available as soon as possible.***  
***Until then, please perform manual installation as described here***

### Download
Use the following command to download the latest release of SARS2seq and move to the newly downloaded `SARS2seq/` directory:

```bash
git clone https://github.com/RIVM-bioinformatics/SARS2seq.git; cd SARS2seq; git checkout tags/$(git tag --sort=committerdate | tail -1) >> /dev/null
```

### Installation
**Before you install SARS2seq, make sure [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) is installed on your system and functioning properly!**

1. Create the required conda-environment and install the necessary dependencies.  
    Copy and paste the code-snippet below in order to create the new conda-environment and directly activate it.  
    
    `conda create --name SARS2seq -c conda-forge mamba python=3.7 -y; conda activate SARS2seq; mamba env update -f mamba-env.yaml`

    You can also use the following snippet if the code-snippet above didn't work for you:  
    `conda env create -f env.yml && conda activate SARS2seq`  
    
    **The "SARS2seq" environment should now be active**  

2. You can now actually install SARS2seq to your system with: `pip install .`

    If that command didn't work for you then use the following: `python setup.py install`

SARS2seq is now installed!  
You can start SARS2seq from anywhere on your system as long as the SARS2seq conda-environment is active.  
You can also use SARS2seq in a different conda-environment as long as the software dependencies match.

### Auto-updating

Since version `v0.3.1`, SARS2seq is able to update itself to the latest released version.  
This makes it easier for everyone to use the latest available version without having to manually check the GitHub releases.

During first-usage of SARS2seq you will be asked to configure your preffered settings for auto-updating. (enabled/disabled/ask-everytime).

If you wish to run SARS2seq without the (auto-)updater checking for new releases, then add the `--skip-updates` flag to your command. In this case you will **not** be notified if there is a new release available, nor will your version of SARS2seq be updated.  

---

## Preparing input primer file

In order to get optimal results, make sure the fasta headers in your fasta file with primers is formatted properly.
Please make sure the fasta headers for your primers are formatted in the following format:

`>{primer-name}_{primer-number}_{orientation}`

Orientation keywords for forward primers are: LEFT/PLUS/POSITIVE  
Orientation keywords for reverse primers are: RIGHT/MINUS/NEGATIVE  

Below is an example of formatted primer names from the [ArticV3](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V3) sequencing protocol:

```Markdown
>nCoV-2019_1_LEFT  
ACCAACCAACTTTCGATCTCTTGT  
>nCoV-2019_1_RIGHT  
CATCTTTAAGATGTTGACGTGCCTC  
>nCoV-2019_2_LEFT  
CTGTTTTACAGGTTCGCGACGT  
>nCoV-2019_2_RIGHT  
TAAGGATCAGTGCCAAGCTCGT
```


If your protocol has alternative primers then make sure the fasta header contains the "alt" keyword in the following format:

`>{primer-name}_{primer-number}_alt_{orientation}`  

Make sure the "alt" keyword is in the middle and not at the end of the fasta header.

Below is an example of formatted primer names from the [ArticV3](https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V3) sequencing protocol with alternative primers included:

```Markdown
>nCoV-2019_13_LEFT  
TCGCACAAATGTCTACTTAGCTGT  
>nCoV-2019_13_RIGHT  
ACCACAGCAGTTAAAACACCCT  
>nCoV-2019_14_LEFT  
CATCCAGATTCTGCCACTCTTGT  
>nCoV-2019_14_alt_LEFT  
TGGCAATCTTCATCCAGATTCTGC  
>nCoV-2019_14_RIGHT  
AGTTTCCACACAGACAGGCATT  
>nCoV-2019_14_alt_RIGHT  
TGCGTGTTTCTTCTGCATGTGC  
>nCoV-2019_15_LEFT  
ACAGTGCTTAAAAAGTGTAAAAGTGCC  
>nCoV-2019_15_alt_LEFT  
AGTGCTTAAAAAGTGTAAAAGTGCCT  
>nCoV-2019_15_RIGHT  
AACAGAAACTGTAGCTGGCACT  
>nCoV-2019_15_alt_RIGHT  
ACTGTAGCTGGCACTTTGAGAGA
```  

Formatting your input primer file as above ensures the best results for your analysis.

---
## Running an analysis

Please see the command line help for a quick explanation of every possible argument: `sars2seq -h`

You can start an analysis with a command such as the following:
```bash
sars2seq \
    --input {path/to/FastQ-files} \
    --output {path/to/desired-output} \
    --primers {path/to/primers.fasta/NONE} \
    --platform {nanopore/illumina/iontorrent} \
    --amplicon-type {end-to-end/end-to-mid} \
    --threads {threads}
``` 

Here, threads refers to the amount of threads that SARS2seq may use on your LOCAL machine. If you're using SARS2seq on a HPC/cluster then these threads will only be used during the pre-processing steps.  
If you're using SARS2seq on your local machine then the given amount of threads will act as a 'ceiling' of threads during analysis.

If your protocol does not use primers then set "NONE" for the primers flag with: `--primers NONE`

The `--amplicon-type` flag is used to clarify the length of of the sequenced read in regards to the amplicon.  
"end-to-end" means that the sequenced read covers the full length of the amplicon, meaning that the primer-sequence is present at both ends of a read.
"end-to-mid" means that the sequenced read *partially* covers the length of the amplicon, meaning that the primer-sequence is present at only one end of a read.  
Please see [this link](https://rivm-bioinformatics.github.io/AmpliGone/latest/amplicon-types/) for a more detailed explanation of the terminology.  

Please check your sequencing and laboratory setup to ensure the best results.

---
## Authors

* Florian Zwagemaker
* Dennis Schmitz
* Karim Hajji
* Annelies Kroneman

## Acknowledgements

* Harry Vennema
* Dirk Eggink
* Jeroen Cremer
* Jeroen Laros 
* Robert Verhagen
* Erwin van Wieringen