import argparse
import sys

import pandas as pd


def GetArgs(sysargs):

    args = argparse.ArgumentParser()

    args.add_argument(
        "-index", type=str, required=True, help="Input nucleotide distribution index"
    )
    args.add_argument(
        "-primers", type=str, required=True, help="Input with primer coordinates"
    )
    args.add_argument("-output", type=str, required=False, help="Output file")

    flags = args.parse_args(sysargs)
    return flags


def split_frames(df):
    left = ["LEFT", "PLUS", "POSITIVE", "FORWARD"]
    right = ["RIGHT", "MINUS", "NEGATIVE", "REVERSE"]

    leftdf = pd.DataFrame(columns=df.columns)
    rightdf = pd.DataFrame(columns=df.columns)

    for x in df.itertuples():
        if any(l in x[1] for l in left) is True:
            leftdf = leftdf.append(
                pd.DataFrame(
                    {"name": x.name, "start": x.start, "stop": x.stop}, index=[0]
                )
            )
        if any(r in x[1] for r in right) is True:
            rightdf = rightdf.append(
                pd.DataFrame(
                    {"name": x.name, "start": x.start, "stop": x.stop}, index=[0]
                )
            )

    leftdf.reset_index(inplace=True)
    rightdf.reset_index(inplace=True)

    leftdf = leftdf.drop(columns=["index"])
    rightdf = rightdf.drop(columns=["index"])

    return leftdf, rightdf


def remove_keyword(prname):
    keywords = [
        "LEFT",
        "RIGHT",
        "PLUS",
        "MINUS",
        "POSITIVE",
        "NEGATIVE",
        "FORWARD",
        "REVERSE",
    ]
    sname = prname.split("_")
    for y, z in enumerate(sname):
        if z in keywords:
            sname.pop(y)
    return "_".join(sname)


def remove_alt_keyword(df):
    for x, y in enumerate(df["name"]):
        if "alt" in y:
            z = y.split("_")
            for a, b in enumerate(z):
                if "alt" in b:
                    z.pop(a)
            z = "_".join(z)
            df.at[x, "name"] = z
    return df


def index_to_remove_stops(one, indexone, two, indextwo):
    stop_one = one.get("stop")
    stop_two = two.get("stop")
    if stop_one > stop_two:
        return indextwo
    if stop_two > stop_one:
        return indexone
    return None


def index_to_remove_starts(one, indexone, two, indextwo):
    start_one = one.get("start")
    start_two = two.get("start")
    if start_one < start_two:
        return indextwo
    if start_two < start_one:
        return indexone
    return None


def remove_alt_primer_l(df):
    xx = df.to_dict(orient="records")
    to_rm = []
    lastindex = list(enumerate(xx))[-1][0]
    for a, x in enumerate(xx):
        if a != lastindex:
            if xx[a].get("name") == xx[a + 1].get("name"):
                rm_indx = index_to_remove_stops(xx[a], a, xx[a + 1], a + 1)
                if rm_indx is not None:
                    to_rm.append(rm_indx)
    filtereddf = df.drop(to_rm)
    return filtereddf


def remove_alt_primer_r(df):
    xx = df.to_dict(orient="records")
    to_rm = []
    lastindex = list(enumerate(xx))[-1][0]
    for a, x in enumerate(xx):
        if a != lastindex:
            if xx[a].get("name") == xx[a + 1].get("name"):
                rm_indx = index_to_remove_starts(xx[a], a, xx[a + 1], a + 1)
                if rm_indx is not None:
                    to_rm.append(rm_indx)
    filtereddf = df.drop(to_rm)
    return filtereddf


def FindUniqueCoords(df):
    dd = df.to_dict(orient="records")
    startingpoint = {}
    endingpoint = {}
    lastindex = list(enumerate(dd))[-1][0]
    firstindex = list(enumerate(dd))[0][0]
    for x, v in enumerate(dd):
        t_end = v.get("rightstart")
        if x != firstindex:
            s = dd[x - 1].get("rightstart")
        else:
            s = v.get("leftstop")
        if x != lastindex:
            end_override = dd[x + 1].get("leftstop")
        else:
            end_override = None
        if end_override is not None:
            if end_override in range(s, t_end):
                primerstart = s
                primerend = end_override
            else:
                primerstart = s
                primerend = t_end
        else:
            primerstart = s
            primerend = t_end
        startingpoint[primerstart] = v.get("name")
        endingpoint[primerend] = v.get("name")

    startdf = (
        pd.DataFrame.from_dict(startingpoint, orient="index")
        .reset_index()
        .rename(columns={0: "name", "index": "unique_start"})
    )
    enddf = (
        pd.DataFrame.from_dict(endingpoint, orient="index")
        .reset_index()
        .rename(columns={0: "name", "index": "unique_end"})
    )
    df = pd.merge(df, startdf, on="name", how="inner")
    df = pd.merge(df, enddf, on="name", how="inner")

    return df


def GetOverlapCoords(df):
    data = df.to_dict(orient="records")
    x = []

    lastindex = list(enumerate(data))[-1][0]
    for k, v in enumerate(data):
        if k != lastindex:
            s = v.get("unique_end")
            e = data[k + 1].get("unique_start")
            x.append([i for i in [*range(s, e)]])

    return [a + 1 for b in x for a in b]


def FilterIndex(ix, coords):
    return ix[ix.index.isin(coords)]


def ReadIndex(f):
    return pd.read_csv(f, index_col=0, compression="gzip")


def ReadPrimers(f):
    return pd.read_csv(f)


def primerframes(input):
    data = ReadPrimers(input)

    if data.empty:
        return None

    l, r = split_frames(data)
    l["name"] = l["name"].apply(remove_keyword)
    r["name"] = r["name"].apply(remove_keyword)

    ll = remove_alt_primer_l(remove_alt_keyword(l))
    rr = remove_alt_primer_r(remove_alt_keyword(r))
    return ll, rr


if __name__ == "__main__":
    flags = GetArgs(sys.argv[1:])

    primers = primerframes(flags.primers)
    indx = ReadIndex(flags.index)
    if primers is None:
        data = pd.DataFrame(columns=list(indx.columns))
    else:
        left, right = primers

        data = FilterIndex(
            indx,
            GetOverlapCoords(
                FindUniqueCoords(
                    pd.merge(left, right, on="name", how="inner")
                    .rename(
                        columns={
                            "start_x": "leftstart",
                            "stop_x": "leftstop",
                            "start_y": "rightstart",
                            "stop_y": "rightstop",
                        }
                    )
                    .drop_duplicates(subset="name")
                )
            ),
        )

    data.to_csv(flags.output, sep=",", index=True, compression="gzip")
