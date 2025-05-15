import pandas as pd

path = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Tuca/CI/"

centromere = "centromere2.bed" 
chrsize = "chrsize"

df = pd.read_csv(path+centromere, sep="\t", header=None)
df["middle"] =  (df[2] + df[1])//2
df2 = pd.read_csv(path+chrsize, sep="\t", header=None)

df_merged = pd.merge(df.loc[:,[0,"middle"]], df2, on=0)
df_merged["half"] = df_merged[1]//2

def middle_arm(middle, half, len):
    if middle > half:
        arm = len - middle
    else:
        arm = middle
    return arm

df_merged["parm"] = df_merged.apply(lambda x: middle_arm(x["middle"], x["half"], x[1]), axis=1)

df_merged["CI"] = df_merged["parm"]/df_merged[1]*100
df_merged["morphology"] = df_merged["CI"].apply(lambda x: "A/T" if x < 25 else "M/SM")
df_merged.sort_values(1, ascending=True, inplace=True)
df_merged.to_csv(path+"CI.tsv", sep="\t", index=False)
