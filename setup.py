import sys

if sys.version_info.major != 3 or sys.version_info.minor < 7:
    print("Error: you must execute setup.py using Python 3.7 or later")
    sys.exit(1)
    
from setuptools import find_packages, setup

from SARS2seq.version import __version__

with open("README.md", "rb") as readme:
    DESCR = readme.read().decode()


setup(
    name="SARS2seq",
    author='Florian Zwagemaker, Dennis Schmitz',
    author_email='rivm-bioinformatics@rivm.nl',
    license='AGPLv3',
    version=__version__,
    packages=find_packages(),
    scripts=[
        'SARS2seq/workflow/workflow.smk',
        'SARS2seq/workflow/directories.py'],
    package_data={'SARS2seq': ['workflow/envs/*', 'workflow/scripts/*', 'workflow/files/*']},
    install_requires=[
        'gnureadline>=8.0.0',
        'biopython>=1.78',
        'snakemake>=6.0.5'
    ],
    entry_points={"console_scripts": [
        'sars2seq = SARS2seq.SARS2seq:main',
        'SARS2seq = SARS2seq.SARS2seq:main']},
    keywords=[],
    include_package_data=True,
    zip_safe=False
)