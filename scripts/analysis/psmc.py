import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.ticker import FuncFormatter
path = '/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/popsize/'
missio = "missioneira_size/"
tuca = "tuca_size/"
files = '.out.msmc2.final.txt'
chrlist = list(range(1,25))


### parametros poblacionales
mu = 3.5e-9
gen = 3
mugen = mu*gen

def plot(chrlist, path, files, chrgray, axis):
    mean = pd.DataFrame()
    for chr in chrlist:
        read = path + "LG" + str(chr)  + files
        df = pd.read_csv(read, sep = "\t")
        x = df["left_time_boundary"] / mugen
        y = (1 / df["lambda"]) / (2 * mu)
        tmpdf =  pd.DataFrame({"x": x, "y": y, "chr":chr}).reset_index()

        if chr == chrgray:
            alpha=1
            color = "gray"
        else:
            alpha = 0.1
            color = "gray"
        sns.lineplot(x=x, y=y, label=str(chr), alpha=alpha, color= color, ax = axis) 
        mean = pd.concat([mean, tmpdf], ignore_index=True)
        
    
    x = mean.groupby("index").x.median()
    y = mean.groupby("index").y.median()
    sns.lineplot(x=x, y=y, label=str(chr), color= "red", ax = axis) 
    axis.set_xscale("log")
    axis.set_ylim(0, 300000)
    axis.set_xlim(1000, 400000)
    axis.set_xlabel("Years ago")
    axis.set_ylabel("Effective population size")
    axis.get_legend().set_visible(False)
    formatter = FuncFormatter(lambda y, _: f'{int(y/1000)}K')
    axis.yaxis.set_major_formatter(formatter)
    return mean



def plotcoal(chrlist, path, files, chrgray, axis):
    mean = pd.DataFrame()
    for chr in chrlist:
        read = path + "LG" + str(chr)  + files
        df = pd.read_csv(read, sep = "\t")
        x = df["left_time_boundary"] / mugen
        y = df["lambda"] 
        tmpdf =  pd.DataFrame({"x": x, "y": y, "chr":chr}).reset_index()

        if chr == chrgray:
            alpha=1
            color = "gray"
        else:
            alpha = 0.1
            color = "gray"
        sns.lineplot(x=x, y=y, label=str(chr), alpha=alpha, color= color, ax = axis) 
        mean = pd.concat([mean, tmpdf], ignore_index=True)
        
    
    x = mean.groupby("index").x.median()
    y = mean.groupby("index").y.median()
    sns.lineplot(x=x, y=y, label=str(chr), color= "red", ax = axis) 
    axis.set_xscale("log")
    axis.set_ylim(0, 7000)
    axis.set_xlim(0, 500000)
    axis.set_xlabel("Years ago")
    axis.set_ylabel("Scaled coalescence rate λ(t)")
    axis.get_legend().set_visible(False)
    formatter = FuncFormatter(lambda y, _: f'{int(y/1000)}K')
    axis.yaxis.set_major_formatter(formatter)
    return mean

mu = 3.5e-9
gen = 3
mugen = mu*gen

fig, axes = plt.subplots(2, 2, figsize=(10, 10))

missioneira = plot(chrlist, path + missio, files,23, axes[0][1])
missioneira = plotcoal(chrlist, path + missio, files,23, axes[0][0])

df = pd.read_csv(path + missio + "LG1.out.msmc2.final.txt", sep = "\t")
plot(chrlist, path + tuca, files, 3, axes[1][1])
plotcoal(chrlist, path + tuca, files, 3, axes[1][0])

# Add row titles
fig.text(0.5, 0.92, '$C. missioneira$', ha='center', fontsize=14)
fig.text(0.5, 0.48, '$C. tuca$', ha='center', fontsize=14)

plt.savefig("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/popsize/msmc2.svg", dpi=300)
# ##### cross coalescent
# # Constants
# mu = 3.5e-9
# gen = 3

# # Load data
# crossPopDat = pd.read_csv("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/popsize/LG3.combined.out", delim_whitespace=True)

# # Calculate time and relative cross-coalescence rate
# time_years = crossPopDat['left_time_boundary'] / mu * gen
# rccr = 2 * crossPopDat['lambda_01'] / (crossPopDat['lambda_00'] + crossPopDat['lambda_11'])

# # Plot
# plt.figure(figsize=(8, 5))
# plt.step(time_years, rccr, where='pre')
# plt.xlim(1000, 500000)
# plt.ylim(0, 1)
# plt.xlabel("Years ago")
# plt.ylabel("Relative cross-coalescence rate")
# plt.xscale("log")  # optional, depending on MSMC2 output
# plt.grid(True)
# plt.show()
