import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
path = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Tuca/repeats/Chr/"

final_kimura = pd.DataFrame()
final_len = pd.DataFrame()

for number in range(1,25):
    file = "Ctuca_"+ str(number) + ".divsum.csv"
    df = pd.read_csv(path + file, sep="\t", names=["Class","Repeat","absLen","wellCharLen","Kimura%"],skiprows=1)
    # df = df[~df["Class"].str.contains("Satellite", case=False)]
    kimura_temp = df.drop(columns=["Class", "absLen","wellCharLen"]).set_index("Repeat").T
    kimura_temp["Chr"] = str(number)
    len_temp = df.drop(columns=["Class", "absLen","Kimura%"]).set_index("Repeat").T
    len_temp["Chr"] = str(number)
    kimura_temp.set_index("Chr", inplace=True)
    len_temp.set_index("Chr", inplace=True) 
    final_kimura = pd.concat([final_kimura,kimura_temp], join="outer", axis=0)
    final_len = pd.concat([final_len,len_temp], join="outer", axis=0)

final_len

above100k = final_len.loc[:, (final_len > 150000).any()]
above100k.drop(columns="combined", inplace=True)
above100k = above100k.fillna(0)
above100k.reset_index(inplace=True)
above100k["Chr"] = above100k.Chr.astype(int)
# order = [17,23,8,9,7,19,15,4,10,13,3,21,16,22,18,20,11,2,14,12,6,24,5,1]
order = list(range(1,25))
above100k['Chr'] = pd.Categorical(above100k['Chr'], categories=order, ordered=True)
sorted_df = above100k.sort_values('Chr').drop(columns="Chr")

final_kimura = final_kimura.loc[:,sorted_df.columns]
final_kimura = final_kimura.astype(float)
final_kimura = final_kimura[(final_kimura<50) & (final_kimura>0)]


import matplotlib as mpl
import numpy as np

cmap = mpl.colormaps.get_cmap('viridis').reversed()
cmap.set_bad("k")

# final_kimura.reset_index(inplace=True)
# final_kimura["Chr"] = final_kimura.Chr.astype(int)
# order = [17,23,8,9,7,19,15,4,10,13,3,21,16,22,18,20,11,2,14,12,6,24,5,1]
# final_kimura['Chr'] = pd.Categorical(final_kimura['Chr'], categories=order, ordered=True)
# final_kimura = final_kimura.sort_values('Chr').drop(columns="Chr")

# fig, ax = plt.subplots(1,2,figsize=(25, 10), sharey=True, sharex=True)
fig, ax = plt.subplots(1,1,figsize=(10, 10))

sns.heatmap(sorted_df, ax=ax,  xticklabels=True)
ax.set_yticklabels(list(range(1,25)))
# sns.heatmap(final_kimura, cmap=cmap, ax=ax[1],  xticklabels=True)

plt.savefig(path+"heatmapchr.svg")

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

path = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Tuca/repeats/Chr/"


fig, ax= plt.subplots(6,2, figsize=(20,30))
axes = ax.flatten()
order = ["RC","Satellite","LTR","DNA","LINE", "Unknown"]
for number in range(1,13):
    file = "Ctuca_"+ str(number) + ".csv"
    df = pd.read_csv(path + file, delim_whitespace=True).set_index("Div")
    df=df.stack().reset_index()
    df = df[df.Div < 30]
    df[0] = df[0] / 1000000
    landscape = sns.histplot(x="Div", weights=0, hue="level_1", data=df, multiple="stack", bins=30, hue_order=order, palette=color_list_7, alpha=1, edgecolor=None, ax=axes[number-1])
    landscape.spines['top'].set_visible(False)
    landscape.spines['right'].set_visible(False)
    landscape.set_xlabel("Divergence from consensus (%)")
    landscape.set_ylabel("Abundance (Mb)")
    landscape.set_title("Chr "+str(number))


# plt.savefig(path + "landscape1.png", dpi=300)


fig, ax= plt.subplots(6,2, figsize=(20,30))
axes = ax.flatten()
order = ["RC","Satellite","LTR","DNA","LINE", "Unknown"]

index = 0
for number in range(13,25):
    file = "Ctuca_"+ str(number) + ".csv"
    df = pd.read_csv(path + file, delim_whitespace=True).set_index("Div")
    df=df.stack().reset_index()
    df = df[df.Div < 30]
    df[0] = df[0] / 1000000
    landscape = sns.histplot(x="Div", weights=0, hue="level_1", data=df, multiple="stack", bins=30, hue_order=order, palette=color_list_7, alpha=1, edgecolor=None, ax=axes[index])
    landscape.spines['top'].set_visible(False)
    landscape.spines['right'].set_visible(False)
    landscape.set_xlabel("Divergence from consensus (%)")
    landscape.set_ylabel("Abundance (Mb)")
    landscape.set_title("Chr "+str(number))
    index += 1


# plt.savefig(path + "landscape2.png", dpi=300)

file = "Ctuca_4.csv"
df = pd.read_csv(path + file, delim_whitespace=True).set_index("Div").T.reset_index()
df = df[df["index"].str.contains("Unknown")]
df.set_index("index", inplace=True)
df = df.iloc[:,0:2]
df = df.loc[df.sum(axis=1) > 10000].T

