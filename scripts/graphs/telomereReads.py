import pandas as pd
import seaborn as sns
import subprocess
import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt


def process_tsv(df):
    df["window"] = df.groupby('id').cumcount() + 1
    df["TTAGGG"] = df["forward_repeat_number"] + df["reverse_repeat_number"]
    return df


def plot_end(df, path, basename, ax=None, save=True):
    if ax == None:
        fig, ax = plt.subplots(1,1)
    sns.lineplot(df, x="window", y="TTAGGG", ax=ax)
    ax.set_ylabel("TTAGGG count per 10Kb")
    ax.set_xlabel("10Kb windows from end")
    ax.set_xlim(1,25)
    ax.set_xticks([0, 5, 10, 15, 20])

    if save:
        plt.savefig(path +"/fig/"+ basename + ".png")
        plt.close()

def blast(reads, reference):
    command = ["makeblastdb", "-in", reads, "-dbtype", "nucl"]
    subprocess.call(command)
    command2 = ["blastn","-query",reference,"-db",reads, '-outfmt', "6 qseqid sseqid pident length slen qstart qend sstart send sstrand","-out", reads + ".out","-dust", "no"]
    subprocess.call(command2)
    return reads + ".out"

def telomere_id(reads, path, basenamse):
    command = ["tidk", "search", "-s", "TTAGGG", "-w", "10000", "-o", basename, "--dir", path, reads]
    subprocess.call(command)
    df = pd.read_csv(path + basename + "_telomeric_repeat_windows.tsv", sep="\t")
    return df

def read_blast(outfile):
    blast = pd.read_csv(outfile, sep="\t", header=None, names=["qseqid","sseqid","pident","length","slen","qstart","qend", "sstart", "send", "strand"])
    blast = blast[blast["length"]>4000]
    blast.sort_values(["length","pident"], inplace=True)
    blast = blast.groupby("sseqid").first().reset_index()
    blast.loc[blast.qstart > blast.qend, "strand"] = "minus"
    blast_dict =  blast.loc[:,["sseqid","strand"]].set_index("sseqid").to_dict()["strand"]
    
    return blast_dict, blast 

def writefasta(blast_dict,reads, end):
    filtered_records = []
    for record in SeqIO.parse(reads, "fasta"): 
        if record.id in blast_dict:
            if blast_dict[record.id] =="minus" or end == "3":
                record.seq = record.seq.reverse_complement()
            filtered_records.append(record)
        
    SeqIO.write(filtered_records, reads + "inverse.fasta", "fasta")
    return reads + "inverse.fasta"


def rightend(x, end):
    if end == "3":
        if x.strand == "plus":
            return x["slen"] - x["send"]
        elif x.strand == "minus":
            return x["sstart"]
    else: 
        if x.strand == "plus":
            return x["sstart"]
        elif x.strand == "minus":
            return x["slen"] - x["sstart"]

path = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Tuca/telomere/FQ/"
end3 = list(range(1,25))
max_dict3 = {}
for chr_end in end3:
    end = "3"
    basename = "LG" + str(chr_end) + "_" + end
    reference = path +"/ref/" + basename
    reads = path + "/fasta/"+ basename + ".fasta"
    blastout = blast(reads, reference)
    blast_dict, blast_df = read_blast(blastout)
    maximo = blast_df.apply(lambda s: rightend(s, end), axis=1).max()
    max_dict3[  str(chr_end) + "_" + end ] = maximo
    reads_ = writefasta(blast_dict, reads, end)
    telomere = telomere_id(reads_, path+"tsv/", basename)
    telomere  = process_tsv(telomere)
    plot_end(telomere, path, basename)

max_dict3



end3 = list(range(1,25))
max_dict5 = {}
for chr_end in end3:
    end = "5"
    basename = "LG" + str(chr_end) + "_" + end
    reference = path + "/ref/" + basename
    reads = path +"/fasta/"+ basename + ".fasta"
    blastout = blast(reads, reference)
    blast_dict, blast_df = read_blast(blastout)
    maximo = blast_df.apply(lambda s: rightend(s, end), axis=1).max()
    max_dict5[  str(chr_end) + "_" + end ] = maximo
    reads_ = writefasta(blast_dict, reads, end)
    telomere = telomere_id(reads_, path+"tsv/", basename)
    telomere  = process_tsv(telomere)
    plot_end(telomere, path, basename)
    

max_dict5




## just plot
########### loop to make a panel
path = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Tuca/telomere/FQ/"
end3 = list(range(1,25))
max_dict5 = {}
fig, axes = plt.subplots(6, 4, sharex=True, constrained_layout=True, figsize=(9,9))
axes = axes.ravel()  
for chr_end in end3:
    end = "5"
    basename = "LG" + str(chr_end) + "_" + end
    reference = path + "/ref/" + basename
    reads = path +"/fasta/"+ basename + ".fasta"
    blastout = reads +".out" 
    blast_dict, blast_df = read_blast(blastout)
    maximo = blast_df.apply(lambda s: rightend(s, end), axis=1).max()
    max_dict5[  str(chr_end) + "_" + end ] = maximo
    reads_ = writefasta(blast_dict, reads, end)
    telomere = pd.read_csv(path +"tsv/"+ basename + "_telomeric_repeat_windows.tsv", sep="\t")
    telomere  = process_tsv(telomere)
    plot_end(telomere, path, basename, axes[chr_end-1])
    axes[chr_end-1].set_ylabel("")
    axes[chr_end-1].set_xlabel("")
    nreads = telomere.id.nunique()
    axes[chr_end-1].set_title("Chr " + str(chr_end)+"\nn="+str(nreads), fontsize=10)
fig.supxlabel("10Kb windows from 5' end", fontsize=14)
fig.supylabel("TTAGGG count per window", fontsize=14)
plt.savefig(path + "/fig/" + basename +".png")

max_dict3 = {}
fig, axes = plt.subplots(6, 4, sharex=True, constrained_layout=True, figsize=(8,8))
axes = axes.ravel()  
for chr_end in end3:
    end = "3"
    basename = "LG" + str(chr_end) + "_" + end
    reference = path + "/ref/" + basename
    reads = path +"/fasta/"+ basename + ".fasta"
    blastout = reads +".out" 
    blast_dict, blast_df = read_blast(blastout)
    maximo = blast_df.apply(lambda s: rightend(s, end), axis=1).max()
    max_dict3[  str(chr_end) + "_" + end ] = maximo
    reads_ = writefasta(blast_dict, reads, end)
    telomere = pd.read_csv(path +"tsv/" + basename + "_telomeric_repeat_windows.tsv", sep="\t")
    telomere  = process_tsv(telomere)
    plot_end(telomere, path, basename, axes[chr_end-1])
    axes[chr_end-1].set_ylabel("")
    axes[chr_end-1].set_xlabel("")
    nreads = telomere.id.nunique()
    axes[chr_end-1].set_title("Chr " + str(chr_end)+"\nn="+str(nreads), fontsize=10)
fig.supxlabel("10Kb windows from 3' end", fontsize=14)
fig.supylabel("TTAGGG count per window", fontsize=14)
plt.savefig(path + "/fig/" + basename +".png")





