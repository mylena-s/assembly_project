#!/usr/bin/env python
import csv
from collections import defaultdict
import pandas as pd

# Initialize counts
counts = defaultdict(lambda: defaultdict(int))

# Input and output file names
input_file = '/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/processed/merged_clean.vcf'
output_file = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/processed/shared.vcf"

# Process the file
with open(input_file, "r") as file:
    for line in file:
        # Skip empty or malformed lines
        if not line.strip():
            continue
        
        # Parse the line
        columns = line.strip().split("\t")
        if len(columns) < 8:  # Ensure the line has enough columns
            continue
        
        # Extract relevant fields
        info_field = columns[7]  # INFO field is in the 8th column (index 7)
        info_data = {key: value for key, value in (item.split("=") for item in info_field.split(";") if "=" in item)}
        
        # Get SVTYPE and SUPP_VEC
        svtype = info_data.get("SVTYPE", "NA")
        supp_vec = info_data.get("SUPP_VEC", "NA")
        # Count occurrences if SUPP_VEC is relevant
        if supp_vec in {"10","01", "11"}:
            counts[svtype][supp_vec] += 1

# Write results to a CSV file
with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    # Write header
    writer.writerow(["SVTYPE", "SUPP_VEC=10","SUPP_VEC=01", "SUPP_VEC=11"])
    for svtype in sorted(counts.keys()):
        writer.writerow([svtype, counts[svtype]["10"], counts[svtype]["01"], counts[svtype]["11"]])


df = pd.read_csv(input_file, sep="\t", comment="#", header=None)
df["type"] = df[2].apply(lambda x: x.split(".")[1])
# exclusive 01, homozygous
df[df[7].str.contains("SUPP_VEC=01")][df[10].str.contains("1/1")].groupby(0).size()

# exclusive 01, heterozygous
df[df[7].str.contains("SUPP_VEC=01")][df[10].str.contains("0/1")].groupby(0).size()

# exclusive 10, heterozygous
df[df[7].str.contains("SUPP_VEC=10")][df[9].str.contains("0/1")].groupby(0).size()
# shared, heterozygous 1 
df[df[7].str.contains("SUPP_VEC=11")][df[9].str.contains("0/1")].groupby(0).size()

# shared, heterozygous 1 and heterozigous in 2
df[df[7].str.contains("SUPP_VEC=11")][df[9].str.contains("0/1")][df[10].str.contains("0/1")].groupby(0).size()
# shared, heterozygous 1 and homozygous in 2
df[df[7].str.contains("SUPP_VEC=11")][df[9].str.contains("0/1")][df[10].str.contains("1/1")].groupby(0).size()

# what trans species polymorphism?
df[df[7].str.contains("SUPP_VEC=11")][df[9].str.contains("0/1")][df[10].str.contains("0/1")].groupby("type").size()
