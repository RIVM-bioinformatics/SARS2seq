import sys
from datetime import date

import pandas as pd

samplename, nx_vs, pn_vs, nx_typ, pn_typ, outfile = sys.argv[1:]

nextclade = pd.read_csv(nx_typ, sep=";").to_dict(orient="list")
pangolin = pd.read_csv(pn_typ, sep=",")

if len(pangolin.index) < 1:
    pangolin = {
        "version": [None],
        "lineage": [None],
        "scorpio_call": [None],
        "qc_status": [None],
    }

else:
    pangolin = pangolin.to_dict(orient="list")

now = date.today().strftime("%d-%m-%Y")

with open(nx_vs, "r") as nf:
    nx_version = nf.read().split(" ")[-1].rstrip()

with open(pn_vs, "r") as pf:
    pn_version = pf.read().split(" ")[-1].rstrip()

with open(outfile, "w") as out:
    out.write(
        f"{samplename}\t{now}\t{pn_version}\t{nx_version}\t{pangolin.get('version')[0]}\t{pangolin.get('lineage')[0]}\t{nextclade.get('clade')[0]}\t{pangolin.get('scorpio_call')[0]}\t{pangolin.get('qc_status')[0]}\t{nextclade.get('qc.overallStatus')[0]}\n"
    )
