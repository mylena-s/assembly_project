import pandas as pd
import seaborn as sns
import statsmodels.api as sm 
import numpy as np
import pybedtools

path = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Oniloticus/centromere"
entropy = pd.read_csv(path + '/Chr14.fasta_entropy.tsv.wig.bw.bigwig.bedgraph', sep="\t", header=None)
complexity = pd.read_csv(path + '/Chr14.fasta_complexity.tsv.wig.bw.bigwig.bedgraph', sep="\t", header=None)

sns.distplot(entropy[3])
sns.distplot(complexity[3])

ethereshold = entropy[3].mean() - 1.5* entropy[3].std()
entropy = entropy[entropy[3] < ethereshold]
entropy.to_csv(path + "/Chr14.entropy_subset.bedgraph", sep="\t", header=False, index=False)

cthereshold = complexity[3].mean() - 2 * complexity[3].std()
complexity = complexity[complexity[3] < cthereshold]
complexity.to_csv(path + "/Chr14.complexity_subset.bedgraph", sep="\t", header=False, index=False)

c = pybedtools.BedTool(path + "/Chr14.complexity_subset.bedgraph")
e = pybedtools.BedTool(path + "/Chr14.entropy_subset.bedgraph")

c = c.merge(d=10000)
e = e.merge(d=10000)

def len_filter(feature, L):
    "Returns True if feature is longer than L"
    return len(feature) > L

intersection = c+e
intersection = intersection.filter(len_filter, L=100000)
intersection.saveas(path + '/Chr14_complexity_entropy_intersection.bed')
