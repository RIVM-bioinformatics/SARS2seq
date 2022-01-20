import sys

from packaging import version as vv
from setuptools import find_packages, setup

from SARS2seq.version import __version__

if sys.version_info.major != 3 or sys.version_info.minor < 7:
    print("Error: you must execute setup.py using Python 3.7 or later")
    sys.exit(1)

try:
    import conda
except SystemError:
    sys.exit(
        """
Error: conda could not be accessed.
Please make sure conda is installed and functioning properly before installing SARS2seq
"""
    )

try:
    import snakemake
except SystemError:
    sys.exit(
        """
Error: SnakeMake could not be accessed.
Please make sure SnakeMake is installed properly before installing SARS2seq
"""
    )

if vv.parse(snakemake.__version__) < vv.parse("6.0"):
    sys.exit(
        f"""
The installed SnakeMake version is older than the minimally required version:

Installed SnakeMake version: {snakemake.__version__}
Required SnakeMake version: 6.0 or later

Please update SnakeMake to a supported version and try again
"""
    )

with open("README.md", "rb") as readme:
    DESCR = readme.read().decode()


setup(
    name="SARS2seq",
    author="Florian Zwagemaker, Dennis Schmitz, Karim Hajji, Annelies Kroneman",
    author_email="ids-bioinformatics@rivm.nl",
    license="AGPLv3",
    version=__version__,
    packages=find_packages(),
    scripts=["SARS2seq/workflow/workflow.smk", "SARS2seq/workflow/directories.py"],
    package_data={
        "SARS2seq": ["workflow/envs/*", "workflow/scripts/*", "workflow/files/*"]
    },
    install_requires=["urllib3>=1.26", "biopython>=1.78", "drmaa==0.7.9", "fpdf==1.7.2", "pyyaml==6.0"],
    entry_points={
        "console_scripts": [
            "sars2seq = SARS2seq.SARS2seq:main",
            "SARS2seq = SARS2seq.SARS2seq:main",
        ]
    },
    keywords=[],
    include_package_data=True,
    zip_safe=False,
)
