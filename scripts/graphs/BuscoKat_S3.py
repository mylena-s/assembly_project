import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import colors

path = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figure1/"
### set global variables for plotting
params_dict= {"axes.titlesize" :24,
"axes.labelsize" : 20,
"lines.linewidth" : 5,
"lines.markersize" : 10,
"xtick.labelsize" : 20,
"ytick.labelsize" : 20,
"legend.frameon":False,
"legend.title_fontsize": 18}

matplotlib.rcParams.update(params_dict)


cmap5 = colors.ListedColormap(["#E9950E","#D2520D","#006CAA","#101D4B","#00A074"])

cmap4 = colors.ListedColormap(["#D2520D","#00A074","#006CAA","#E9950E"])

### 1. process statistic table
## Contiguity
stats = pd.read_csv(path + "CBFH00195.stats", sep="\t", header=None, names=["Statistic", "value"]).set_index("Statistic").T
stats = stats.T.drop_duplicates().T
stat_cols = ["# scaffolds", "Total scaffold length", "Scaffold N50", "Scaffold L50", "Scaffold N90","Scaffold L90","GC content %"]
round_cols = ["Total scaffold length", "Scaffold N50", "Scaffold N90"]
stats = stats[stat_cols].astype(float)
stats.loc[:,round_cols] = round(stats[round_cols] / 1000000, 1)
stats.loc[:,"GC content %"] = round(stats["GC content %"], 1)
stats[["Scaffold L90", "Scaffold L50"]] = stats.loc[:,["Scaffold L90", "Scaffold L50"]].apply(np.int64)
chromosomes = 24
stats["Scaffold to chromosome ratio"] = stats["# scaffolds"].apply(lambda x: round(x / chromosomes, 2))
stats["Scaffold N50/L50"] = stats.apply(lambda x: str(x["Scaffold N50"]) + "/" + str(int(x["Scaffold L50"])),axis=1)
stats["Scaffold N90/L90"] = stats.apply(lambda x: str(x["Scaffold N90"]) + "/" + str(x["Scaffold L90"]),axis=1)
stats = stats.drop(["Scaffold N50", "Scaffold L50", "Scaffold N90", "Scaffold L90"], axis=1).T
stats.rename(columns={"value":"Statistics"}, inplace=True)
### 2. Completeness

fig, ax = plt.subplots(1,2, figsize=(20,4), gridspec_kw={'width_ratios': [5,4], "hspace":2})
# kmer completeness
kmer = pd.read_csv(path+"kat-main.mx", skiprows=12, sep="\s", header=None)
xlim = 80
kmer = kmer.iloc[0:xlim, 0:5]
kmer_plot = kmer.plot.area(stacked=True, ylim=(0,35000000), ylabel="Distinct k-mers",linewidth=0, xlabel="k-mer multiplicity", ax=ax[0],cmap=cmap5)
kmer_plot.spines['top'].set_visible(False)
kmer_plot.spines['right'].set_visible(False)
kmer_plot.set_xticks(range(0,xlim, 10))
completeness = round(pd.read_json(path + "kat.dist_analysis.json")["completeness"]["k"],2)
kmer_plot.set_title("k-mer completeness")
kmer_plot.legend(title="Times in assembly",prop={'size': 15}, loc="lower center", frameon=False,  bbox_to_anchor=(0.5, -0.75), ncols=3)

# gene completeness
busco_genome = pd.read_json(path + "assembly_short_summary.specific.actinopterygii_odb10.01_busco_out_CBFH00195.json")
statistics = ["Complete", "Single copy", "Multi copy", "Fragmented", "Missing","n_markers"]
odb_v = busco_genome.loc["name", "lineage_dataset"] 
busco_genome = busco_genome.loc[statistics, "results"]


#the same for gene annotations braker
busco_annotations_braker = pd.read_json(path + "/Annotations_LongestIso/short_summary.specific.actinopterygii_odb10.braker_longestiso_busco.json")
busco_annotations_braker = busco_annotations_braker.loc[statistics, "results"]

#the same for gene annotations galba
busco_annotations_galba = pd.read_json(path + "/Annotations_LongestIso/short_summary.specific.actinopterygii_odb10.galba_longestiso_busco.json")
busco_annotations_galba = busco_annotations_galba.loc[statistics, "results"]


#the same for gene annotations galba
busco_annotations_toga = pd.read_json(path + "/Annotations_LongestIso/short_summary.specific.actinopterygii_odb10.toga_longestiso_busco.json")
busco_annotations_toga = busco_annotations_toga.loc[statistics, "results"]


busco = pd.DataFrame()
busco["Genome"] = busco_genome 
busco["Braker"] = busco_annotations_braker # (Change df)
busco["Galba"] = busco_annotations_galba # (Change df)
busco["TOGA"] = busco_annotations_toga # (Change df)

busco = busco.T
busco = busco.loc[["TOGA", "Galba", "Braker", "Genome"]] 
busco_plot = busco.iloc[:, 1:-1].plot.barh(stacked=True,title ="Gene completeness", ylabel="",linewidth=0, cmap=cmap4,xlabel="Percentage of genes (n="+str(busco_genome.n_markers)+")", xlim=(0,100),ax=ax[1], width=0.8)
busco_plot.spines['top'].set_visible(False)
busco_plot.spines['left'].set_visible(False)
busco_plot.spines['right'].set_visible(False)
busco_plot.set_xticks(range(0,101, 20))
ax[1].yaxis.tick_right()
plt.legend(title=odb_v, loc="lower center",prop={'size': 15}, ncol=2, frameon=False, bbox_to_anchor=(0.5, -0.75))

plt.subplots_adjust(hspace=5)

plt.savefig(path + "output/completeness_iso.svg", dpi=300)

### 3. Correctness
inspector = pd.read_csv(path + "summary_statistics", sep ="\t")
insp_cols = ['Mapping rate /%', 'Split-read rate /%','Depth',
       'Structural error', 'Expansion', 'Collapse', 'Haplotype switch',
       'Inversion', 'Small-scale assembly error /per Mbp', 'QV']

inspector=inspector.loc[insp_cols,:].round(2)
inspector.rename(columns={"Statics of contigs:" :"Statistics"}, inplace=True)
pd.concat([stats,inspector]).to_csv(path+"output/processed_statistics.csv", sep="\t")
