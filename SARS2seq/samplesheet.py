"""
Write the samplesheets
"""

import re
import yaml
import os

def illumina_sheet(inputdir, sheet):
    illuminapattern = re.compile("(.*)(_|\.)R?(1|2)(?:_.*\.|\..*\.|\.)f(ast)?q(\.gz)?")
    samples = {}
    for dirname, subdir, filename in os.walk(inputdir):
        for files in filename:
            fullpath = os.path.join(dirname, files)
            match = illuminapattern.fullmatch(files)
            if match:
                sample = samples.setdefault(match.group(1), {})
                sample["R{}".format(match.group(3))] = str(fullpath)
    with open(sheet, "w") as samplesheet:
        yaml.dump(samples, samplesheet, default_flow_style=False)
    samplesheet.close()
    pass

def nanopore_sheet(inputdir, sheet):
    nanoporepattern = re.compile("(.*)\.f(ast)?q(\.gz)?")
    samples = {}
    for dirname, subdir, filename in os.walk(inputdir):
        for files in filename:
            fullpath = os.path.join(dirname, files)
            match = nanoporepattern.fullmatch(files)
            if match:
                samples.setdefault(match.group(1), fullpath)
    with open(sheet, "w") as samplesheet:
        yaml.dump(samples, samplesheet, default_flow_style=False)
    samplesheet.close()
    pass


def WriteSampleSheet(input, platform):
    if platform == 'illumina':
        illumina_sheet(input, "samplesheet.yaml")
    if platform == 'nanopore':
        nanopore_sheet(input, "samplesheet.yaml")
    
    samplesheet = os.getcwd() + "/samplesheet.yaml"
    return samplesheet