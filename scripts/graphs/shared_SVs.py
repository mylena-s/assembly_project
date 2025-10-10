#!/usr/bin/env python
import csv
from collections import defaultdict
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
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
df = df[~df[9].str.contains("1/1")]
df["SUPP_VEC"] = df[7].str.extract(r"SUPP_VEC=(\d{2})")
df["$C. missioneira$ genotype"] = df[10].str.extract(r"^(\d\/\d)")
df["$C. tuca$ genotype"] = df[9].str.extract(r"^(\d\/\d)")
df=df.fillna(0)
df["SVLEN"] = df[7].str.extract(r"SVLEN=([-+]?\d+)")
df["SVLEN"] = df["SVLEN"].apply(lambda x: np.abs(int(x)))
df2 = df.copy()
df = df[(df["type"]=="INV") | (df["type"]=="DUP")]
df = df[df[0]!= "LG10"]
all_sp_df = pd.DataFrame(np.nan,index=["Count","Lenght"], columns=["$C. tuca$","$C. missioneira$","$A. citrinellus$","$O. niloticus$"])
all_sp_df.loc["Count","$C. tuca$"] = df[(df["SUPP_VEC"]=="10") |(df["SUPP_VEC"]=="11")].count().unique()[0]
all_sp_df.loc["Lenght","$C. tuca$"] = df[(df["SUPP_VEC"]=="10") |(df["SUPP_VEC"]=="11")]["SVLEN"].sum()
all_sp_df.loc["Count","$C. missioneira$"] = df[(df["SUPP_VEC"]=="01") |(df["SUPP_VEC"]=="11")].count().unique()[0]
all_sp_df.loc["Lenght","$C. missioneira$"] = df[(df["SUPP_VEC"]=="01") |(df["SUPP_VEC"]=="11")]["SVLEN"].sum()

all_sp_df.loc["Count","$A. citrinellus$"] = 256 + 10266 
all_sp_df.loc["Lenght","$A. citrinellus$"] = 142260256 + 18987282 
all_sp_df.loc["Count","$O. niloticus$"] = 161 + 6941 
all_sp_df.loc["Lenght","$O. niloticus$"] = 149958659 + 9686288 
amp_onil = 153380894 + 11828557

syntenic = 505875690 - 415064392
cmap = sns.color_palette("Blues", as_cmap=True)

fig, axes = plt.subplots(1,2, figsize=(8,4))
sns.heatmap(pd.DataFrame(all_sp_df.loc["Lenght",]), annot=True, annot_kws={"size": 12},fmt=",.0f",cmap=cmap, cbar=False, ax=axes[0])

pivot = df2.groupby(["$C. tuca$ genotype", "$C. missioneira$ genotype"])["SVLEN"].sum().unstack(fill_value=0)

# pivot = df2.groupby(["$C. tuca$ genotype", "$C. missioneira$ genotype"]).size().unstack(fill_value=0)
pivot=pivot.replace(0, np.nan)
pivot=pivot.rename(columns={0:"0/0"},index={0:"0/0"})
pivot = pivot.T

custom_annotations = pd.DataFrame([
    ["-", "Intra-\nspecies", "Inter-\nspecies"],
    ["Intra-\nspecies", "Trans-\nspecies", "Trans/intra-\nspecies"]]).T


annot = np.empty_like(pivot, dtype=object)
for i in range(pivot.shape[0]):
    for j in range(pivot.shape[1]):
        annot[i, j] = f"{pivot.iloc[i, j]:,.0f}\n\n{custom_annotations.iloc[i, j]}"
pivot.iloc[0,1]

sns.heatmap(pivot, annot=annot, annot_kws={"size": 11},fmt="",cmap=cmap, ax=axes[1], cbar=False)
axes[1].tick_params(axis="x")
for _, spine in axes[1].spines.items():
        spine.set_visible(True)
for _, spine in axes[0].spines.items():
        spine.set_visible(True)

axes[0].set_xlabel("[Inversions + Duplications]")

axes[0].set_xticklabels(["$C. tuca$"])
axes[1].set_xticklabels(["Hom.","Het."])
axes[1].set_yticklabels(["Hom. ref.","Het.","Hom."])

axes[0].set_title("Structural divergence\n in Cichlidae (bp)", fontsize=14)
axes[1].set_title("Structural diversity\n in $Crenicichla$ (bp)", fontsize=14)
plt.tight_layout()

plt.savefig("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figure_SVacrosstime/FigDivergence.svg", dpi=500)




#### estimate lengths distributions from syri vcf

comp1=pd.read_csv("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figure2/asm20_10k/CBFH00195_synteny_Amphilopus_syntenysyri.vcf", sep="\t", comment="#", header=None)
comp2=pd.read_csv("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figure2/asm20_10k/Amphilopus_synteny_Oniloticus_syntenysyri.vcf", sep="\t", comment="#", header=None)
comp1 = comp1[comp1[4].isin(["<DUP>","<INVDP>", "<INV>","<TRANS>","<INVTR>"])]
comp2 = comp2[comp2[4].isin(["<DUP>","<INVDP>", "<INV>","<TRANS>","<INVTR>"])]
comp1[4] =  comp1[4].str.replace("<INVDP>","<DUP>")
comp1[4] =  comp1[4].str.replace("<INVTR>","<TRANS>")
comp2[4] =  comp2[4].str.replace("<INVDP>","<DUP>")
comp2[4] =  comp2[4].str.replace("<INVTR>","<TRANS>")
comp1.groupby(4).describe().to_csv("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figure2/asm20_10k/tuca_midas.description")
comp2.groupby(4).describe().to_csv("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figure2/asm20_10k/midas_tilapia.description")