from Bio import SeqIO
import re

def ContainsSpecials(seq):
    
    chars = re.compile("[^actgumrwsykvhdbnACTGUMRWSYKVHDBN-]")
    
    if (chars.search(seq)== None):
        return False
    else:
        return True
    
def IsValidFasta(inputfile):
    results = []
    for record in SeqIO.parse(inputfile, "fasta"):
        results.append(ContainsSpecials(str(record.seq)))
    
    if any(results) is True:
        return False
    else:
        return True