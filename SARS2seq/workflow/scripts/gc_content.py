import argparse
import pandas as pd
from Bio import SeqIO
from Bio import SeqUtils
from decimal import *

pd.set_option("display.max_rows", None, "display.max_columns", None)


# conda install -c bioconda seqiolib
# didn't work right away

# conda install biopython 
# after this it did work but complained about env being unstable.

# python ../gc_content.py --uniqlocs Test_NotTranspose_R826_B96_4802157154.csv --key R826_B96_4802157154 --sequence R826_B96_4802157154_cov_ge_100.fa

# newline
# python /mnt/scratch_dir/landmanf/coverage_plot_new/testdir/gc_content.py --uniqlocs uniqlocs_newoutput_R826_B96_4802157154.csv --key R826_B96_4802157154 --sequence R826_B96_4802157154_cov_ge_100.fa --entire_table entire_table_newoutput_R826_B96_4802157154.csv --primers primers.fasta


args = argparse.ArgumentParser()

args.add_argument(
    "--uniqlocs",
    metavar="File",
    type=str,
    help="Output file with unique start and end location",
    required=True,
)

args.add_argument(
    "--key",
    metavar="String",
    help="Sample ID",
    type=str,
    required=True,
)

args.add_argument(
    "--sequence",
    metavar="File",
    help="Consensus sequence fasta file to determine QC content",
    type=str,
    required=True,
)

args.add_argument(
    "--entire_table",
    metavar="File",
    type=str,
    help="Output entire table so I have primer end and start positions that got removed",
    required=True,
)

args.add_argument(
    "--primers",
    metavar="File",
    type=str,
    help="input file with primers as given by AmpliGone",
    required=True,
)

flags = args.parse_args()

record = SeqIO.read(flags.sequence, "fasta") # It is zero indexed
record_list = list(SeqIO.parse(flags.sequence, "fasta"))
getcontext().prec = 4 # Sets the decimal for GC content value

# Make a list of ranges, make a new list to fill up the ranges in between then with that set of ranges determine the sequences and give those to seqio and write to dataframe
# loop through range of entire sequence, get the unique start and end then get ranges based on those known values.

uniqlocs = pd.read_csv(
    flags.uniqlocs, sep=",", names=["index", "name", "unique_start", "unique_end"]
)
uniqlocs_vanilla = uniqlocs.drop(columns="index")
uniqlocs = uniqlocs.drop(columns="index")
uniqlocs["gc_content"] = ""

entire_table = pd.read_csv(
    flags.entire_table, sep=",", names=["index", "name", "leftstart", "leftstop", "rightstart", "rightstop","unique_start", "unique_end"]
)


# What I really would like to do is directly create a new column named uniq_sequence by calling record[start:end] but can't somehow.
# Some futile attempts:
# uniqlocs["record2"] = record[(uniqlocs["unique_start"]):(uniqlocs["unique_end"])]
# uniqlocs["record2"] = record[(uniqlocs["unique_start"].astype(int)):(uniqlocs["unique_end"].astype(int))]
# uniqlocs["uniq_test"] = [f"{x}:{y}" for x, y in zip(uniqlocs["unique_start"], uniqlocs["unique_end"])]
# print(record[(40):(2010)].seq)
# print(SeqUtils.GC(record[340:450].seq))
# print(record[340:450].seq) # This prints the sequence and that's what you want



## CALCULATE THE GC CONTENT ON THE INITIAL DATAFRAME WITH UNIQUE COORDS

# Could also just make one with the key names and gc_content and then merge??
# Could make this a function depending on the dataframe that it receives? So that it can work for dataframe with only unique amplicon and the dataframe with all.
for x in range(len(uniqlocs)):
    # Think I still have to correct for 0 index here
    start = int(uniqlocs.loc[x, ("unique_start")])
    end = int(uniqlocs.loc[x, ("unique_end")])
    sequence = record[(start):(end)].seq
    gc_content = SeqUtils.GC(sequence)
    uniqlocs.loc[x, ("gc_content")] = round(gc_content, 2) # Only prints 2 decimals if there's no zero already

# print(uniqlocs)
# print(record[54:342].seq) # This prints the sequence and that's what you want


# CALCULATE THE NON UNIQUE PARTS
final_index = (len(uniqlocs_vanilla))-1
# Lengte is 98, x gaat tot 97
for x in range(len(uniqlocs_vanilla)):
    name = uniqlocs_vanilla.loc[x, ("name")]
    if x == 0:
        pre_non_uniq_start = 1
        pre_non_uniq_end = int(uniqlocs_vanilla.loc[x, ("unique_start")]) -1
        uniqlocs_vanilla = uniqlocs_vanilla.append({'name': f"pre_{name}", "unique_start":pre_non_uniq_start, "unique_end":pre_non_uniq_end}, ignore_index=True)
        continue

    if x == final_index:
        pre_non_uniq_start = int(uniqlocs_vanilla.loc[x-1, ("unique_end")]) +1
        pre_non_uniq_end = int(uniqlocs_vanilla.loc[x, ("unique_start")]) -1
        uniqlocs_vanilla = uniqlocs_vanilla.append({'name': f"pre_{name}", "unique_start":pre_non_uniq_start, "unique_end":pre_non_uniq_end}, ignore_index=True)

        post_non_uniq_start = int(uniqlocs_vanilla.loc[x, ("unique_end")]) +1
        post_non_uniq_end = len(record)
        uniqlocs_vanilla = uniqlocs_vanilla.append({'name': f"post_{name}", "unique_start":post_non_uniq_start, "unique_end":post_non_uniq_end}, ignore_index=True)

    else:
        pre_non_uniq_start = int(uniqlocs_vanilla.loc[x-1, ("unique_end")]) +1
        pre_non_uniq_end = int(uniqlocs_vanilla.loc[x, ("unique_start")]) -1
        uniqlocs_vanilla = uniqlocs_vanilla.append({'name': f"pre_{name}", "unique_start":pre_non_uniq_start, "unique_end":pre_non_uniq_end}, ignore_index=True)


## THIS NOW HAS THE NON UNIQUE PARTS ADDED IN THE DATAFRAME 
uniqlocs_vanilla["gc_content"] = ""
for x in range(len(uniqlocs_vanilla)):
    # Think I still have to correct for 0 index here
    start = int(uniqlocs_vanilla.loc[x, ("unique_start")])
    end = int(uniqlocs_vanilla.loc[x, ("unique_end")])
    sequence = record[(start):(end)].seq
    gc_content = SeqUtils.GC(sequence)
    uniqlocs_vanilla.loc[x, ("gc_content")] = round(gc_content, 2)


uniqlocs_vanilla_sorted = uniqlocs_vanilla.sort_values(by="unique_start")
uniqlocs_vanilla_sorted = uniqlocs_vanilla_sorted.reset_index(drop=True)
# print(uniqlocs_vanilla_sorted)

# Make an primer dictionary to look for the primer sequence when later looping through the entire table
def ConvertPrimerstoDict(primerfile):
    primers_dictionary = {}
    with open(primerfile,'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines):
            if line.startswith('>'):
                key = line.strip('>').strip()
                sequence = lines[i + 1].strip()
                primers_dictionary[key] = sequence
    return(primers_dictionary)

def listToString(s): 
    str1 = "" 
    for ele in s: 
        str1 += ele  
    return str1 

def Countstring(strng, char):
    count = 0
    for x in strng:
        if x == (str(char)):
            count = count + 1
    return(count)

def sequence_compare(seq_a, seq_b):
        len1= len(seq_a)
        len2= len(seq_b)
        mismatches = []
        for pos in range (0,min(len1,len2)) :
            if seq_a[pos] != seq_b[pos]:
                mismatches.append('|')
            else:
                mismatches.append(' ')
        return(listToString(mismatches))

primer_mismatch_list = []
# Loop through the entire table with start/end of left and right primers and compare its sequence to actual primer sequence.
for x in range(len(entire_table)):
    primers_dictionary = ConvertPrimerstoDict(flags.primers)
    left_start = int(entire_table.loc[x, ("leftstart")])
    left_end = int(entire_table.loc[x, ("leftstop")])
    right_start = int(entire_table.loc[x, ("rightstart")])
    right_end = int(entire_table.loc[x, ("rightstop")])
    name = entire_table.loc[x, ("name")]
    sequence_primer_left = str(record[(left_start):(left_end)].seq)
    sequence_primer_right = str((record[(right_start):(right_end)].seq).reverse_complement())
    dictsearchkey = name.split("_")[0] + "-" + name.split("_")[1] + "_" + str(int(name.split("_")[2]))

# For the LEFT primers:
    if (dictsearchkey+'_LEFT') in primers_dictionary.keys():
        seq_l = primers_dictionary[(dictsearchkey+'_LEFT')]
        if (dictsearchkey+'_alt_LEFT') in primers_dictionary.keys():
            alt_seq_l = primers_dictionary[(dictsearchkey+'_alt_LEFT')]

        # If neither of the primer matches 100%:
            if (Countstring(sequence_compare(sequence_primer_left, alt_seq_l), '|')) > 0 and (Countstring(sequence_compare(sequence_primer_left, seq_l), '|')) > 0:
                value_list = (Countstring(sequence_compare(sequence_primer_left, alt_seq_l), '|')),(Countstring(sequence_compare(sequence_primer_left, seq_l), '|'))
                minimum_value = min(value_list)
                min_index = value_list.index(minimum_value)
                if min_index == 0:
                    primer_mismatch_list.append([f"{name}_alt_LEFT", minimum_value, len(sequence_primer_left), sequence_primer_left, alt_seq_l])
                else:
                    primer_mismatch_list.append([f"{name}_LEFT", minimum_value, len(sequence_primer_left), sequence_primer_left, seq_l])


        # If one of the primer matches 100%:
            else:
                value_list = (Countstring(sequence_compare(sequence_primer_left, alt_seq_l), '|')),(Countstring(sequence_compare(sequence_primer_left, seq_l), '|'))
                minimum_value = min(value_list)
                min_index = value_list.index(minimum_value)
                if min_index == 0:
                    primer_mismatch_list.append([f"{name}_alt_LEFT", minimum_value, len(sequence_primer_left), sequence_primer_left, alt_seq_l])
                else:
                    primer_mismatch_list.append([f"{name}_LEFT", minimum_value, len(sequence_primer_left), sequence_primer_left, seq_l])

    # No alt primer for this primer:
        else:
            primer_mismatch_list.append([f"{name}_LEFT", (Countstring(sequence_compare(sequence_primer_left, seq_l), '|')), len(sequence_primer_left), sequence_primer_left, seq_l])

# For the RIGHT primers:
    if (dictsearchkey+'_RIGHT') in primers_dictionary.keys():
        seq_r = primers_dictionary[(dictsearchkey+'_RIGHT')]
        if (dictsearchkey+'_alt_RIGHT') in primers_dictionary.keys():
            alt_seq_r = primers_dictionary[(dictsearchkey+'_alt_RIGHT')]

        # If neither of the primer matches 100%:
            if (Countstring(sequence_compare(sequence_primer_right, alt_seq_r), '|')) > 0 and (Countstring(sequence_compare(sequence_primer_right, seq_r), '|')) > 0:
                value_list = (Countstring(sequence_compare(sequence_primer_right, alt_seq_r), '|')),(Countstring(sequence_compare(sequence_primer_right, seq_r), '|'))
                minimum_value = min(value_list)
                min_index = value_list.index(minimum_value)

                if min_index == 0:
                    primer_mismatch_list.append([f"{name}_alt_RIGHT", minimum_value, len(sequence_primer_right), sequence_primer_right, alt_seq_r])

                else:
                    primer_mismatch_list.append([f"{name}_RIGHT", minimum_value, len(sequence_primer_right), sequence_primer_right, seq_r])

        # If one of the primer matches 100% determine which one and append appropiate sequence:
            else:
                value_list = (Countstring(sequence_compare(sequence_primer_right, alt_seq_r), '|')),(Countstring(sequence_compare(sequence_primer_right, seq_r), '|'))
                minimum_value = min(value_list)
                min_index = value_list.index(minimum_value)
                if min_index == 0:
                    primer_mismatch_list.append([f"{name}_alt_RIGHT", minimum_value, len(sequence_primer_right), sequence_primer_right, alt_seq_r])
                else:
                    primer_mismatch_list.append([f"{name}_RIGHT", minimum_value, len(sequence_primer_right), sequence_primer_right, seq_r])

        else:
            primer_mismatch_list.append([f"{name}_RIGHT", (Countstring(sequence_compare(sequence_primer_right, seq_r), '|')), len(sequence_primer_right), sequence_primer_right, seq_r])

primer_mismatch_df = pd.DataFrame(primer_mismatch_list, columns = ['name', 'mismatch', 'length', 'isolate_primer_seq', 'reference_primer_seq'])

primer_mismatch_list.to_csv('primer_mismatch_df.csv', sep=",", index=True, header=False)

# https://pandas.pydata.org/pandas-docs/stable/user_guide/style.html#Export-to-Excel
