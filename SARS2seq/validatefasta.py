import re

from Bio import SeqIO


def ContainsSpecials(seq):

    chars = re.compile("[^actgumrwsykvhdbnACTGUMRWSYKVHDBN-]")

    if chars.search(seq) == None:
        return False
    else:
        return True


def IsValidFasta(inputfile):
    if inputfile == "NONE":
        return True
    results = []
    for record in SeqIO.parse(inputfile, "fasta"):
        results.append(ContainsSpecials(str(record.seq)))

    if any(results) is True:
        return False
    else:
        return True
