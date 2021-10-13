# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## 0.2.2
### added
- This changelog
### changed
#### **bugfixes**
- Fixed issue [#10](https://github.com/RIVM-bioinformatics/SARS2seq/issues/10)  
#### **other**
- Improved `setup.py` to ensure SARS2seq runs in an environment with a supported SnakeMake version.  
    Checks SnakeMake version and exits with proper message if the installed SnakeMake version will cause problems.
- Pinned `libffi` versions in workflow-environments to version `3.3` as we found that later versions of `libffi` may cause problems during environment-setup on some HPC-systems.
- The HTML snakemake report will now be created without command-line output.
- The SnakeMake will now try to complete a failed rule 3 more times before it is actually marked as a failed rule. This is to compensate for filesystem latency on some HPC-systems.
