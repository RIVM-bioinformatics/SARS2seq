import argparse
from datetime import date

import pandas as pd

# pd.set_option('display.max_rows', 500)

arg = argparse.ArgumentParser()

arg.add_argument("--key", metavar="<key>", type=str, required=False)

arg.add_argument(
    "--typing_aggs",
    type=str,
    nargs="+",
    default="1",
    metavar="File1 File2 File3",
    required=False,
)

arg.add_argument(
    "--boc",
    type=str,
    metavar="File",
    required=False,
)

arg.add_argument(
    "--coverages",
    type=int,
    nargs="+",
    default=1,
    metavar="1 5 10",
    help="Coverage thresholds to use",
    required=False,
)

arg.add_argument("--output", type=str, metavar="File", required=False)

args = arg.parse_args()

flags = arg.parse_args()


def readboc(f, cols):
    return pd.read_csv(f, sep="\t", header=None, names=cols)


def choose_best_match(df):
    a = df.to_dict(orient="records")
    ind = None
    for x, v in enumerate(a):
        if v.get("boc") >= 97:
            ind = x

    if ind is None:
        for x, v in enumerate(a):
            if v.get("boc") >= 95:
                ind = x

    if ind is None:
        b = {}
        for x, v in enumerate(a):
            b[x] = v.get("boc")
        highest = max(b, key=b.get)
        ind = highest

    return df.loc[[ind]]


def read_match(frame):
    a = frame.to_dict(orient="records")
    f = a[0].get("corresponding_file")

    nf = pd.read_csv(f, header=None, sep="\t")
    nf = nf.rename(
        columns={
            0: "sample",
            1: "date",
            2: "pangoversion",
            3: "nextversion",
            4: "pangodate",
            5: "pangolineage",
            6: "clade",
            7: "scorpio",
            8: "pangoqc",
            9: "nextqc",
        }
    )
    return nf


def getcov(frame):
    a = frame.to_dict(orient="records")
    x = a[0].get("covs")
    return x


if __name__ == "__main__":

    now = date.today().strftime("%d-%m-%Y")

    bocdf = readboc(flags.boc, flags.coverages)

    bocdf = (
        bocdf.transpose()
        .reset_index()
        .rename(columns={"index": "covs", flags.key: "boc"})
    )

    bocdf["corresponding_file"] = flags.typing_aggs

    match = choose_best_match(bocdf)

    best_typing_agg = read_match(match)

    best_typing_agg["coverage_level"] = getcov(match)

    df = best_typing_agg[
        [
            "sample",
            "coverage_level",
            "date",
            "pangoversion",
            "nextversion",
            "pangodate",
            "pangolineage",
            "clade",
            "scorpio",
            "pangoqc",
            "nextqc",
        ]
    ]

    df.to_csv(flags.output, sep="\t", header=None, index=None)
