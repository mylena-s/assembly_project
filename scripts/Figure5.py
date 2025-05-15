import pandas as pd
import seaborn as sns
import csv
import matplotlib.pyplot as plt
from pybedtools import BedTool


# color_list_7 = ['#101D4B', '#006CAA', '#5FB2E6', '#00A074', "#D2520D",'#E9950E', '#FAEC26']
color_list_7 = ["#00A074","#5FB2E6","#006CAA","#101D4B","#861657","#E9950E","#D2520D"]
color_list_8 = color_list_7.copy()
color_list_6 = color_list_7.copy()
color_list_6.pop(1)
color_list_8.reverse()
color_list_8.insert(0, '#d3d3d3')

def analyze_repeats_in_svs(variants_file, repeats_file, repeat_class_file, genome_repeat_bases, genome_total_bases, repeat_bases_list, repeat_types_list):
    # Load BedTool objects
    variants = BedTool(variants_file)
    repeats = BedTool(repeats_file)

    # Intersect repeats with SVs
    columns2 = ['chrom', 'software', 'repeatname', 'start', 'end', 'xx', 'strand', 'id', 'chr2', 'start2', 'end2', 'svid','geno', 'overlap_len']
    repeat_svs = repeats.intersect(variants, wo=True, f=0.5).to_dataframe(names=columns2)

    # Total repeats in SVs
    repeat_bases_in_svs = repeat_svs.overlap_len.sum()

    # Load repeat classification
    repeat_class = pd.read_csv(repeat_class_file, sep="\t")
    repeat_svs["element"] = repeat_svs.id.str.extract(r'Target=([^\s;]+)')
    repeat_svs = pd.merge(repeat_svs, repeat_class, left_on='element', right_on='element', how='left')
    repeat_svs["classification"] = repeat_svs["classification"].fillna("Simple_repeat")

    # Summarize by repeat type
    repeat_lengths = {rtype: repeat_svs.loc[repeat_svs.classification.str.contains(rtype, na=False), "overlap_len"].sum() for rtype in repeat_types}

    # Return results

    return {
        "Repeat type lengths (SVs)": repeat_lengths,
        "Repeat bases in SVs": repeat_bases_in_svs,
        "df" : repeat_svs,
    }


def divsum_splitter(file):    
    file_modif = file.split('.')[0]
    with open(file) as file, open( file_modif + '.csv', 'w') as modif:    
        orig_table = csv.reader(file)
        for row in orig_table:
            for index in row:
                if 'Div' in index:
                    file_csv = csv.writer(modif)
                    file_csv.writerow(row)
                    for row in orig_table:
                        csv.writer(modif)
                        file_csv.writerow(row)
    return file_modif + '.csv'


path = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figure5/"
divsum_splitter(path + 'CBFH00195.divsum')
df = pd.read_csv(path + 'CBFH00195.csv', delim_whitespace=True).set_index("Div")
df=df.stack().reset_index()
df["Class"] = df.level_1.apply(lambda x: x.split("/")[0])
dfClass = df.groupby(["Class", "Div"]).sum().reset_index()
dfClass[1] = dfClass[0] // 1000
dfClass[0] = dfClass[0] / 1000000
# dfClass = dfClass[~dfClass["Class"].str.contains("Unknown")]
dfClass = dfClass[dfClass.Div < 30]
#Do not include simple repeats because count = 0 and unclassified because its divergence does not make sense
order = ["RC","Satellite","LTR","DNA","LINE", "Unknown"]
counts = pd.read_fwf(path+"types_counts.csv",  sep="\s+", header=None).set_index(0)
counts.loc["Non repetitive",4] = 100-counts.sum()[4]
counts.loc["Simple & Low compl"] = counts.loc["Simple repeats"] + counts.loc["Low complexity"]
counts.columns = ["Copy_number","Length","bp","Proportion"]

order2 = ['Non repetitive',
 'Unclassified',
 'LINEs',
 'DNA',
 'LTR',
 'Satellites',
 'Simple & Low compl',
 'Rolling-circles']

counts = pd.DataFrame(counts["Proportion"]).T
counts.rename(columns={"Unknown":"Unclassified"}, inplace=True)
counts=counts.loc[:,order2]

def reindex_df(df, weight_col):
    """expand the dataframe to prepare for resampling
    result is 1 row per count per sample"""
    df = df.reindex(df.index.repeat(df[weight_col]))
    df.reset_index(drop=True, inplace=True)
    return(df)

DF = reindex_df(dfClass, weight_col = 1)
fig, ax= plt.subplots(1,3, figsize=(10, 5), gridspec_kw={'width_ratios': [1,2,1], "wspace":0.5})
barplot = counts.plot.bar(stacked=True, color=color_list_8, ax=ax[0])
barplot.spines['top'].set_visible(False)
barplot.spines['right'].set_visible(False)
barplot.spines['bottom'].set_visible(False)
barplot.set_ylabel("Genome proportion (%)")
plt.tick_params(left = False)
barplot.xaxis.set_visible(False)
ticks = list(range(0,101,20))
barplot.set_yticks(ticks)
sns.move_legend(
    ax[0], "lower center",
    bbox_to_anchor=(2.8, -0.3), ncol=4, title=None, frameon=False,
)

violinplot = sns.violinplot(y='Class', x='Div', data=DF, orient="h", palette=color_list_6, order=order, ax=ax[1], alpha=1, saturation=1, cut=0)
violinplot.spines['top'].set_visible(False)
violinplot.spines['right'].set_visible(False)
violinplot.spines['left'].set_visible(False)
plt.tick_params(left = False)
violinplot.set_ylabel(" ")
violinplot.set_xlabel("Divergence from consensus (%)")

############# SVS

path = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figure_SVs/"
variants = path + 'processed/CBFH00195_sniffles_clean.bed'
genes = BedTool(path + 'data/gene_toga.gtf')  
introns = BedTool(path + 'data/introns_toga.gtf')  
exons = BedTool(path + 'data/exons_toga.gtf')  
repeats = BedTool(path + 'data/repeats_synteny.gff')
repeat_class = path + "data/repeat_classif.csv"
repeat_types = ["LINE", "Satellite", "LTR", "DNA", "Unknown", "RC|Rolling", "simple_repeat|Simple_repeat"]
repeat_bases = [80656088, 21381757, 15259963, 65026066, 113377172, 975454, 17148607]

Ctuca_results3 = analyze_repeats_in_svs(variants, 
                                  repeats, 
                                  repeat_class, 
                                  313825107, 
                                  839694952,
                                  repeat_bases,
                                  repeat_types)

Ctuca_results3["df"]["level1"] = Ctuca_results3["df"].classification.apply(lambda x: x.split("/")[0])
counts = Ctuca_results3["df"].groupby("level1").overlap_len.sum()/ 9521034
counts.loc["Non repetitive"] = 1 - counts.sum()
counts= counts *100
counts = pd.DataFrame(counts).T
counts.rename(columns={"Unknown":"Unclassified","Simple_repeat":"Simple & Low compl", "RC":'Rolling-circles',"LINE":"LINEs","Satellite":"Satellites"}, inplace=True)
counts=counts.loc[:,order2]


barplot = counts.plot.bar(stacked=True, color=color_list_8, ax=ax[2])
barplot.spines['top'].set_visible(False)
barplot.spines['right'].set_visible(False)
barplot.spines['bottom'].set_visible(False)
barplot.set_ylabel("Proportion of SV bases (%)")
barplot.xaxis.set_visible(False)
barplot.get_legend().set_visible(False)
barplot.set_yticks(ticks)

element = Ctuca_results3["df"].groupby("element").overlap_len.sum() / 9521034
element.sort_values(ascending=False).head(10)
df = pd.read_csv(path + "processed/kimura.csv", sep="\t")
list_e = ["rnd-1_family-40","rnd-1_family-26","rnd-1_family-144","rnd-1_family-1","circ10-2070","rnd-1_family-0","rnd-1_family-96","rnd-4_family-1857","rnd-1_family-73","rnd-1_family-230"]
points = df[df["Repeat"].isin(list_e)].loc[:,["Class","Kimura%"]]
points["Class"] = points["Class"].apply(lambda x: x.split("/")[0])

#### points upon violinplot 
sns.scatterplot(y=points['Class'], x=points['Kimura%'].astype(float), color='red', zorder=10, ax=ax[1])
path = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figure5/"
plt.savefig(path + "repeat_plot_alternativa.png", bbox_inches='tight', dpi=300)
