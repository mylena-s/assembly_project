from pybedtools import BedTool
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.proportion import proportions_ztest
import numpy as np
 
path = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figure_SVs/"
variants = path + 'processed/CBFH00195_sniffles_clean.bed'
genes = path + 'data/genes.gtf'  
annotations_file = path + 'data/eggnog_output.emapper.annotations'

# def analyze_variants(variants_file, genes_file):
#     """
#     """

#     # Load BedTool objects
#     variants = BedTool(variants_file)
#     genes = BedTool(genes_file)

#     # Genic variants
#     columns = ['chrom', 'start', 'end', 'name', 'chr2', 'start2', 'end2', 'overlap_len']
#     genic_df = variants.intersect(genes, wo=True).to_dataframe(names=columns).reset_index()

#     return genic_df


# Ctuca_results = analyze_variants(variants, genes)
# Cmissioneira_results = analyze_variants(path + 'processed/CBFH00251_sniffles_clean.bed', genes)
# annotations = pd.read_csv(annotations_file, sep="\t", header=None, skiprows=5, skipfooter=3)

# annotations["gene"] = annotations[0].apply(lambda x: x.split(".")[0])
# df[0].to_csv(path +"Ctuca.list", index=False, header=False)

# df = pd.merge(Ctuca_results, annotations, how="left", left_on="end2", right_on="gene")
# df.groupby("level_3").agg({1:"unique", 2:"unique"}).reset_index().to_csv(path + "gene_annotSVs_tuca.csv", sep="\t")
# df["len"] = df["gene"].dropna().apply(lambda x: len(x))



# df = pd.merge(Cmissioneira_results, annotations, how="left", left_on="end2", right_on="gene")
# df = df.groupby("level_3").agg({"gene":"unique", "chrom":"unique",1:"unique",2:"unique"})
# df["len"] = df["gene"].dropna().apply(lambda x: len(x))


# get gene annotations from toga 
path = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figure_SVs/"
variants = path + 'processed/CBFH00195_sniffles_clean.bed'
genes = path + 'data/transcript_toga.gtf'  

variants = BedTool(variants)
genes = BedTool(genes)
columns = ['chrom', 'start', 'end', 'name', 'genotype','chr2', 'softw','feature','start2', 'end2', 'dot', 'strand', 'info', 'ovlen']
genic_df = variants.intersect(genes, wo=True).to_dataframe(names=columns)
array = genic_df["info"].str.extract(r'transcript_id "([^"]+)"')[0].str.split('.').str[0].unique()
genic_df["gene"] = genic_df["info"].str.extract(r'transcript_id "([^"]+)"')[0].str.split('.').str[0]
genic_df["gene_id"] = genic_df["info"].str.extract(r'gene_id "([^"]+)"')[0]
pd.DataFrame(array).to_csv(path + "tuca_enrichment.txt", header=False, index=False)
genic_df.groupby("end").gene_id.nunique().sort_values()

genic_df[genic_df.end == "Sniffles2.DEL.4D6SA"].gene.unique()
