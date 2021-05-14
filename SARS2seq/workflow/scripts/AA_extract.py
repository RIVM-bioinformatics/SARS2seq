import sys

import gffpandas.gffpandas as gffpd
from Bio import SeqIO
from Bio.Seq import Seq

nucleotides, features, samplename, cov, outdir = sys.argv[1:]

## read nucleotides
seq = str(SeqIO.read(nucleotides, "fasta").seq)

## read features
gffdict = gffpd.read_gff3(features).df.to_dict(orient="index")

for k in gffdict.keys():
    orfname = str(gffdict[k].get("attributes").split(";")[1].split("=")[-1])
    start = int(gffdict[k].get("start"))
    end = int(gffdict[k].get("end"))

    sliced_seq = seq[start - 1 : end].replace("-", "")
    AminoAcids = Seq(sliced_seq).translate(to_stop=True)
    
    with open(f"{outdir}/{orfname}/{samplename}_{cov}.fa", "w") as out:
        out.write(f">{samplename}_ORF-{orfname}_cov_{cov}\n" f"{AminoAcids}\n")
