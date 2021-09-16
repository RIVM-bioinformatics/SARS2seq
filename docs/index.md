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
git clone https://github.com/RIVM-bioinformatics/SARS2seq.git; cd SARS2seq
```

### Installation
**Before you install SARS2seq, make sure [Conda](https://docs.conda.io/projects/conda/en/latest/index.html) is installed on your system and functioning properly!**

1. Create the required conda-environment and install the necessary dependencies.  
    Copy and paste the code-snippet below in order to create the new conda-environment and directly activate it.  
    
    `conda create --name SARS2seq -c conda-forge mamba python=3.7; conda activate SARS2seq; mamba env update -f mamba-env.yaml`

    You can also use the following snippet if the code-snippet above didn't work for you:  
    `conda env create -f env.yml && conda activate SARS2seq`  
    
    **The "SARS2seq" environment should now be active**  

2. You can now actually install SARS2seq to your system with: `pip install .`

    If that command didn't work for you then use the following: `python setup.py install`