import pandas as pd
import seaborn as sns
import csv
import matplotlib.pyplot as plt


color_list_7 = ['#101D4B', '#006CAA', '#5FB2E6', '#00A074', "#D2520D",'#E9950E', '#FAEC26']
color_list_8 = color_list_7.copy()
color_list_6 = color_list_7.copy()
color_list_6.pop(1)
color_list_8.reverse()
color_list_8.insert(0, '#d3d3d3')

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


path = "/home/mylequis/Documents/servidor_bk/home/fhenning/assembly_project/CBFH00195/results/annotation/repeat/04_masking/"
divsum_splitter(path + 'CBFH00195.divsum')

df = pd.read_csv(path + 'CBFH00195.csv', delim_whitespace=True).set_index("Div")
df=df.stack().reset_index()
df["Class"] = df.level_1.apply(lambda x: x.split("/")[0])
dfClass = df.groupby(["Class", "Div"]).sum().reset_index()
dfClass[1] = dfClass[0] // 1000
dfClass[0] = dfClass[0] / 1000000
dfClass = dfClass[~dfClass["Class"].str.contains("Unknown")]
dfClass = dfClass[dfClass.Div < 30]
#Do not include simple repeats because count = 0 and unclassified because its divergence does not make sense
order = ["RC","Satellite","LTR","DNA","LINE"]

counts = pd.read_fwf(path+"types_counts.csv",  sep="\s+", header=None).set_index(0)
counts.loc["Non repetitive"] = counts.sum()
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

fig, ax= plt.subplots(2,1, figsize=(4, 6.5), gridspec_kw={'height_ratios': [3, 1]})
landscape = sns.histplot(x="Div", weights=0, hue="Class", data=dfClass, multiple="stack", bins=30, hue_order=order, palette=color_list_6, alpha=1, edgecolor=None, ax=ax[0])
landscape.spines['top'].set_visible(False)
landscape.spines['right'].set_visible(False)
landscape.set_xlabel("Divergence from consensus (%)")
landscape.set_ylabel("Abundance (Mb)")
ax[0].legend_.remove()
barplot = counts.plot.barh(stacked=True, color=color_list_8, ax=ax[1])
barplot.spines['top'].set_visible(False)
barplot.spines['right'].set_visible(False)
barplot.spines['left'].set_visible(False)
barplot.set_xlabel("Genome proportion (%)")
barplot.yaxis.set_visible(False)
plt.legend(loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.5, -2))
plt.tight_layout()

plt.savefig(path + "repeat_plot.svg", bbox_inches='tight', dpi=300)

############## alternative 
def reindex_df(df, weight_col):
    """expand the dataframe to prepare for resampling
    result is 1 row per count per sample"""
    df = df.reindex(df.index.repeat(df[weight_col]))
    df.reset_index(drop=True, inplace=True)
    return(df)

DF = reindex_df(dfClass, weight_col = 1)
fig, ax= plt.subplots(1,2, figsize=(8, 5), gridspec_kw={'width_ratios': [1,2]})

barplot = counts.plot.bar(stacked=True, color=color_list_8, ax=ax[0])
barplot.spines['top'].set_visible(False)
barplot.spines['right'].set_visible(False)
barplot.spines['bottom'].set_visible(False)
barplot.set_ylabel("Genome proportion (%)")
plt.tick_params(left = False)
barplot.xaxis.set_visible(False)
sns.move_legend(
    ax[0], "lower center",
    bbox_to_anchor=(1.6, -0.3), ncol=4, title=None, frameon=False,
)

violinplot = sns.violinplot(y='Class', x='Div', data=DF, orient="h", palette=color_list_6, order=order, ax=ax[1], alpha=1, saturation=1, cut=0)
violinplot.spines['top'].set_visible(False)
violinplot.spines['right'].set_visible(False)
violinplot.spines['left'].set_visible(False)
plt.tick_params(left = False)
plt.ylabel(" ")
violinplot.set_xlabel("Divergence from consensus (%)")
plt.savefig(path + "repeat_plot_alternativa.svg", bbox_inches='tight', dpi=300)