import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import colors


Headers1 = ["Single","Duplicated","Unexpected","Expected", "Missing","remove"]
Headers2 = ["Consistent","Partial mapping","Fragmented", "Inconsistent","Partial_I", "Fragmented_I", "Contamination","Partial_cont","Fragmented_cont", "Unknown", "remove"]

df_completeness = pd.DataFrame(columns= Headers1)
df_consistency = pd.DataFrame()

def process_line(line, headers, sample, method):
    line=line.replace("[","")
    line=line.replace("\n","")
    line=line.replace("]","")
    df =pd.DataFrame([line.split("%")],columns=headers)
    df = df.apply(lambda x: x.str.replace(',', ''))
    df = df.apply(lambda x: x.str.replace(r'[A-Z]\:', '', regex=True))
    df["Sample"] = sample
    df["Method"] = method
    return df

def readlines(file, sample, method):
    with open(file) as fp:
        for i, line in enumerate(fp):
            if i == 5:
                df1=process_line(line, Headers1, sample, method)
            elif i == 10:
                df2=process_line(line, Headers2, sample, method)
    return df1, df2


def loop(file, sample, method, df_completeness, df_consistency):
    df1, df2 = readlines(file, sample, method)
    df_completeness = pd.concat([df_completeness, df1], ignore_index=True)
    df_consistency = pd.concat([df_consistency, df2], ignore_index=True)
    return df_completeness, df_consistency


sampleDict = {"Ctuca_braker":['/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/FigureS4/CtucaBraker.sum', "C. tuca", "Braker2"],
              "C_tuca_galba":['/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/FigureS4/CtucaGalba.sum', "C. tuca", "GALBA"],
              "C_tuca_toga":['/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/FigureS4/CtucaTOGA.sum', "C. tuca", "TOGA"],
              "A_citrinellus_braker":['/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/FigureS4/AcitrinellusBraker.sum', "A. citrinellus", "Braker2"]}

for key, value in sampleDict.items():
    df_completeness, df_consistency = loop(value[0],  value[1], value[2], df_completeness, df_consistency)





#### plot
params_dict= {"axes.titlesize" :18,
"axes.labelsize" : 16,
"lines.linewidth" : 5,
"xtick.labelsize" : 12,
"ytick.labelsize" : 12,
"legend.frameon":False,
"legend.title_fontsize": 18,
"legend.fontsize": 16}
matplotlib.rcParams.update(params_dict)


cmap4 = colors.ListedColormap(["#40B8E9","#40B8E9","#40B8E9","#8B5FF3","#8B5FF3","#8B5FF3","#E6A001","#384551"])
cmap3 = colors.ListedColormap(["#48BFA2","#CFEE8F","#ED1E5E"])



df_completeness[["Single","Duplicated","Missing"]] = df_completeness[["Single","Duplicated","Missing"]].astype(float)
df_completeness["Sample"] = "$" + df_completeness["Sample"] + "$"+"\n" + df_completeness["Method"]

df_consistency[["Consistent","Partial mapping","Fragmented","Inconsistent","Partial_I", "Fragmented_I","Contamination","Unknown"]] = df_consistency[["Consistent","Partial mapping","Fragmented","Inconsistent","Partial_I", "Fragmented_I","Contamination","Unknown"]].astype(float)
df_consistency["Consistent"] = df_consistency["Consistent"] - df_consistency["Partial mapping"] - df_consistency["Fragmented"] 
df_consistency["Inconsistent"] = df_consistency["Inconsistent"] - df_consistency["Partial_I"] - df_consistency["Fragmented_I"] 
df_consistency["Sample"] = "$" +df_consistency["Sample"] +"$" +"\n" + df_consistency["Method"]

fig, ax = plt.subplots(2,1, figsize=(8,12), gridspec_kw={"hspace":0}, sharex=True)
Completeness=df_completeness.plot.bar(stacked=True,x="Sample",y=["Single","Duplicated","Missing"],rot=90,  width = 0.95, ax=ax[0], cmap=cmap3,xlabel="", ylabel="Conserved HOGs (%)")
Completeness.spines['top'].set_visible(False)
Completeness.spines['right'].set_visible(False)
Completeness.legend(title="Completeness",prop={'size': 12}, loc="upper right", frameon=False,  bbox_to_anchor=(1.27, 1), ncols=1)

Consistency=df_consistency.plot.bar(stacked=True,x="Sample", y=["Consistent","Partial mapping","Fragmented","Inconsistent","Partial_I", "Fragmented_I","Contamination","Unknown"], cmap=cmap4, width = 0.95, ax=ax[1], xlabel="", ylabel="Percentage of proteome (%)")
Consistency.spines['top'].set_visible(False)
Consistency.spines['right'].set_visible(False)
Consistency.invert_yaxis()
patterns = ["","","","",
            "//" , "//","//","//",
            "\\\\","\\\\" ,"\\\\" ,"\\\\" ,
            "" ,"", "" , "",
            "//" , "//","//","//",
            "\\\\","\\\\" ,"\\\\" ,"\\\\" ,
            "" ,"", "" , "",
            "" ,"", "" , ""]

bars = Consistency.patches
for bar, hatch in zip(bars, patterns):
    bar.set_hatch(hatch)

Consistency.legend(title="Taxonomic\nconsistency",
                    prop={'size': 12}, loc="upper right", 
                    fontsize="large", frameon=False,  
                    bbox_to_anchor=(1.27, 1), ncols=1,
                    labels=["Consistent","Partial mapping","Fragmented","Inconsistent","Partial mapping","Fragmented","Contamination","Unknown"])

#### Annotate bars with N genes
geneCount = [33982,26243, 22459, 34688]

Completeness

for i in range(len(geneCount)):
    Completeness.annotate(
            "n = "+str(geneCount[i]),                      # Use `label` as label
            (i, 78),         # Place label at end of the bar
            textcoords="offset points", # Interpret `xytext` as offset in points
            ha='center',
            fontsize=12,               # Horizontally center label
            )                      # Vertically align label differently for
                                        # positive and negative values.


plt.xticks(rotation=0, ha='center')

plt.savefig("/home/mylequis/Documents/Doctorado/Papers/Genome_report/02_FIGURES/FigureS4.svg", bbox_inches='tight')
