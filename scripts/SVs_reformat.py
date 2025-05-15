from pybedtools import BedTool
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.proportion import proportions_ztest
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
from matplotlib.gridspec import GridSpec
import matplotlib.patches as mpatches
import numpy as np
from PIL import Image
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

path = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/"
variants = path + 'processed/CBFH00195_sniffles_clean.bed'
variants2 = path + "processed/CBFH00251_sniffles_clean.bed"


genes = BedTool(path + 'tracks/toga.gene.gtf')  
introns = BedTool(path + 'tracks/toga.introns.gtf')  
exons = BedTool(path + 'tracks/toga.exons.gtf')  
repeats = BedTool(path + 'tracks/repeats.final.gff')
repeat_class = path + "tracks/repeat_classif.csv"


def analyze_variants(variants_file, genes_file, total_genome_bases):
    """
    Analyzes structural variants enrichment in intergenic regions.

    Parameters:
    variants_file (str): Path to the structural variants file (BED format).
    genes_file (str): Path to the genes file (GTF/BED format).
    total_genome_bases (int): Total bases in the genome.

    Returns:
    dict: Contains Z-statistic, p-value, and additional calculations.
    """

    # Load BedTool objects
    variants = BedTool(variants_file)
    genes = BedTool(genes_file)
    genes = BedTool(genes_file)
    genes = BedTool.sort(genes)
    genes = genes.merge()

    # Gene length calculation
    genes_df = genes.to_dataframe()
    genes_df["len"] = genes_df["end"] - genes_df["start"]
    genic_bases = genes_df["len"].sum()

    # Intergenic variants
    intergenic_svs = variants.subtract(genes).to_dataframe()
    intergenic_svs["len"] = intergenic_svs["end"] - intergenic_svs["start"]
    intergenic_bases_in_svs = intergenic_svs["len"].sum()

    # Genic variants
    columns = ['chrom', 'start', 'end', 'name', 'chr2', 'start2', 'end2', 'overlap_len']
    genic_df = variants.intersect(genes, wo=True).to_dataframe(names=columns)
    genic_bases_in_svs = genic_df["overlap_len"].sum()

    # Intergenic bases in the genome
    intergenic_bases = total_genome_bases - genic_bases

    # Total bases in SVs
    variants_df = variants.to_dataframe()
    variants_df["len"] = variants_df["end"] - variants_df["start"]
    total_svs_bases = variants_df["len"].sum()

    # Perform z-test
    count = [intergenic_bases_in_svs, intergenic_bases]
    nobs = [total_svs_bases, total_genome_bases]
    stat, pval = proportions_ztest(count, nobs, alternative='larger')

    # Return results
    return {
        "Z-statistic": stat,
        "p-value": pval,
        "Intergenic bases in SVs": intergenic_bases_in_svs,
        "Genic bases in SVs": genic_bases_in_svs,
        "Total SV bases": total_svs_bases,
        "Total genome bases": total_genome_bases,
        "Intergenic bases (genome)": intergenic_bases,
        "Genic bases (genome)": genic_bases,
        "genes_df": genes_df,
        "genic_df": genic_df
}

def analyze_intronic_exonic_svs(variants_file, introns_file, exons_file, genic_bases_in_svs, genic_bases,genes_file):
    """
    Analyzes structural variants enrichment in intronic and exonic regions.

    Parameters:
    variants_file (str): Path to the structural variants file (BED format).
    introns_file (str): Path to the introns file (GTF/BED format).
    exons_file (str): Path to the exons file (GTF/BED format).
    genic_bases_in_svs (int): Total genic bases affected by SVs.
    genic_bases (int): Total genic bases in the genome.

    Returns:
    dict: Contains results for intronic and exonic enrichment analysis.
    """

    from statsmodels.stats.proportion import proportions_ztest
    from pybedtools import BedTool

    # Load BedTool objects
    variants = BedTool(variants_file)
    exons = BedTool(exons_file)

    # as TOGA does not include introns, if added manually with agat, the script breaks
    # so in case of running with galba or braker, uncomment this and comment the other alternative
    # introns = BedTool(introns_file)
    # introns = BedTool.sort(introns)
    # introns = introns.merge()

    exons = BedTool.sort(exons)
    exons = exons.merge()
    
    genes = BedTool(genes_file)
    genes = BedTool.sort(genes)
    genes = genes.merge()
    
    introns = genes.subtract(exons)
    # Intronic SVs
    columns = ['chrom', 'start', 'end', 'name', 'chr2', 'start2', 'end2', 'overlap_len']
    intronic_svs = variants.intersect(introns, wo=True).to_dataframe(names=columns)
    introns_df = introns.to_dataframe()
    introns_df["len"] = introns_df["end"] - introns_df["start"]
    intronic_len = introns_df["len"].sum()
    intronic_len_sv = intronic_svs["overlap_len"].sum()

    # Exonic SVs
    exonic_svs = variants.intersect(exons, wo=True).to_dataframe(names=columns)
    exons_df = exons.to_dataframe()
    exons_df["len"] = exons_df["end"] - exons_df["start"]
    exonic_len = exons_df["len"].sum()
    exonic_len_sv = exonic_svs["overlap_len"].sum()

    # Intronic z-test
    count_introns = [intronic_len_sv, intronic_len]
    nobs_introns = [genic_bases_in_svs, genic_bases]
    stat_introns, pval_introns = proportions_ztest(count_introns, nobs_introns, alternative='larger')

    # Exonic z-test
    count_exons = [exonic_len_sv, exonic_len]
    nobs_exons = [genic_bases_in_svs, genic_bases]
    stat_exons, pval_exons = proportions_ztest(count_exons, nobs_exons, alternative='smaller')

    # Return results
    return {
        "Intronic Z-statistic": stat_introns,
        "Intronic p-value": pval_introns,
        "Intronic bases in SVs": intronic_len_sv,
        "Intronic bases (genome)": intronic_len,
        "Exonic Z-statistic": stat_exons,
        "Exonic p-value": pval_exons,
        "Exonic bases in SVs": exonic_len_sv,
        "Exonic bases (genome)": exonic_len,
        "exons df": exons_df,
        "intron_svs": intronic_svs
    }



def analyze_repeats_in_svs(variants_file, repeats_file, repeat_class_file, genome_repeat_bases, genome_total_bases, repeat_bases_list, repeat_types_list):
    """
    Analyzes the enrichment of repeats in structural variants (SVs) and the proportions of repeat types.

    Parameters:
    variants_file (str): Path to the variants file (BED format).
    repeats_file (str): Path to the repeats file (GTF/BED format).
    repeat_class_file (str): Path to the repeat classification file (CSV format).
    genome_repeat_bases (int): Total repeat bases in the genome.
    genome_total_bases (int): Total bases in the genome.

    Returns:
    dict: Contains results for repeats enrichment and repeat type proportions.
    """
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
    # repeat_svs["element"] = repeat_svs.id
    repeat_svs = pd.merge(repeat_svs, repeat_class, left_on='element', right_on='element', how='left')
    repeat_svs["classification"] = repeat_svs["classification"].fillna("Simple_repeat")

    # Summarize by repeat type
    repeat_lengths = {rtype: repeat_svs.loc[repeat_svs.classification.str.contains(rtype, na=False), "overlap_len"].sum() for rtype in repeat_types}

    # Observations for Z-test: enrichment of repeats in SVs
    total_svs_bases = variants.to_dataframe()["end"].sub(variants.to_dataframe()["start"]).sum()
    count = [repeat_bases_in_svs, genome_repeat_bases]
    nobs = [total_svs_bases, genome_total_bases]
    stat, pval = proportions_ztest(count, nobs, alternative='larger')

    # Chi-square test for repeat type proportions
    observed = [
         [repeat_lengths[rtype] for rtype in repeat_types_list],  # Counts for SVs
         repeat_bases_list  # Counts for genome
     ]
    chi2, chi_pval, dof, expected = stats.chi2_contingency(observed)

    # Individual Z-tests for each repeat type
    repeat_pvals = []
    for i, rtype in enumerate(repeat_types):
        count = [observed[0][i], observed[1][i]]
        stat, rtype_pval = proportions_ztest(count, nobs, alternative='larger')
        repeat_pvals.append(rtype_pval)

    # Adjust p-values for multiple testing
    adjusted_pvals = multipletests(repeat_pvals, method='fdr_bh')[1]

    # Return results

    return {
        "Repeat type lengths (SVs)": repeat_lengths,
        "Repeat bases in SVs": repeat_bases_in_svs,
        "Repeat enrichment Z-statistic": stat,
        "Repeat enrichment p-value": pval,
        "Chi-square statistic": chi2,
        "Chi-square p-value": chi_pval,
        "Adjusted p-values for repeat types": adjusted_pvals,
        "df" : repeat_svs,
    }





repeat_bases_in_genome = 313825107
total_genome_bases = 839694952

Ctuca_results = analyze_variants(
    variants,
    genes,
    total_genome_bases
)

Ctuca_results2 = analyze_intronic_exonic_svs(
    variants,
    introns, 
    exons,
    Ctuca_results['Genic bases in SVs'], 
    Ctuca_results['Genic bases (genome)'],genes)
Ctuca_results2


repeat_types = ["LINE", "Satellite", "LTR", "DNA", "Unknown", "RC|Rolling", "simple_repeat|Simple_repeat"]
repeat_bases = [80656088, 21381757, 15259963, 65026066, 113377172, 975454, 17148607]

Ctuca_results3 = analyze_repeats_in_svs(variants, 
                                  repeats, 
                                  repeat_class, 
                                  313825107, 
                                  839694952,
                                  repeat_bases,
                                  repeat_types)


C_missioneira_results = analyze_variants(
    variants2,
    genes,
    total_genome_bases
)

C_missioneira_results2 = analyze_intronic_exonic_svs(
    variants2,
    introns,
    exons,
    C_missioneira_results['Genic bases in SVs'],
    C_missioneira_results['Genic bases (genome)'], genes)

C_missioneira_results3 = analyze_repeats_in_svs(
    variants2,
    repeats, 
    repeat_class, 
    313825107, 
    839694952,
    repeat_bases,
    repeat_types)


################
#### Graphs ####
################

def process_variants(file_paths, samples):
    """
    Processes variant files and prepares a dataframe for plotting.

    Parameters:
    file_paths (list): List of file paths to variant BED files.
    samples (list): List of sample names corresponding to the files.

    Returns:
    pd.DataFrame: A concatenated dataframe with processed variant data.
    """
    if len(file_paths) != len(samples):
        raise ValueError("The number of file paths must match the number of samples.")

    dataframes = []
    for file_path, sample in zip(file_paths, samples):
        temp = pd.read_csv(file_path, sep="\t", names=["chr", "start", "end", "name", "genotype"])
        temp["sample"] = sample
        dataframes.append(temp)

    # Concatenate all dataframes
    variants_plot = pd.concat(dataframes)

    # Extract variant type and simplify genotype
    variants_plot["type"] = variants_plot.name.str.extract(r'(INS|DEL|DUP|INV|BND)')
    variants_plot["genotype"] = variants_plot["genotype"].str.split(':').str[0]
    variants_plot["length"] = variants_plot.end - variants_plot.start

    bins = [0, 5000, 10000,50000, float('inf')]
    labels = ['<5 Kbp', '5-10 Kbp', '10-50 Kbp', '50 Kbp+']
    # Recode the lenght variable
    variants_plot['length_category'] = pd.cut(variants_plot['length'], bins=bins, labels=labels, right=False)
    # combine DEL and INS into INDEL
    variants_plot["type"] = variants_plot["type"].str.replace(r'DEL', "INDEL")
    variants_plot["type"] = variants_plot["type"].str.replace(r'INS', "INDEL")


    return variants_plot

def format_axes(axes):
    for i, ax in enumerate(axes):
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.tick_params(axis='x', rotation=0, labelsize=12)
        ax.tick_params(axis='y', labelsize=12)
        ax.set_xlabel("")
        # ax.set_ylabel("")
        ax.legend(frameon=False, fancybox=False)

def plot_expected_vs_observed(ax, title, featureprop, nonfeatureprop, bar_width=0.8):
    categories = ["$C. missioneira$","$C. tuca$", "Exp."]
    # Stacked bar plot
    ax.bar(categories, featureprop, label="Feature", color=["#0070bd","#D2520D","darkgray"], width=bar_width)
    ax.bar(categories, nonfeatureprop, bottom=featureprop, label="Non-Feature", color=["#0070bd","#D2520D","darkgray"], hatch="//", edgecolor='white', width=bar_width)
    ax.set_ylabel(title, size=16)
    ax.set_yticks([])
    ax.set_xticks([])

def plot_stacked_prop(df, sample, ax):
    pivoted = df.melt(id_vars=['sample', 'length_category'], var_name='type', value_name='count')
    pivoted['percentage'] = pivoted.groupby(['sample', 'length_category'])['count'].transform(lambda x: x / x.sum() * 100)
    subset = pivoted[pivoted['sample'] == sample]
    subset = subset[subset['type']!="BND"]
    subset_pivot = subset.pivot(index='length_category', columns='type', values='percentage')
    subset_pivot.sort_values("length_category", ascending=False).plot(kind='barh', stacked=True, ax=ax, width = 0.9, label=sample, color=["#00A074","#5FB2E6", "#D2520D"])

def plot_proportions(ax, title, dataset):
    """
    Calculates proportions and plots observed vs expected proportions for a given category.

    Parameters:
    ax (matplotlib.axes.Axes): Matplotlib axis object to plot on.
    title (str): Title for the plot.
    dataset (dict): Dictionary containing the dataset with counts and totals.
    category (str): Category name to extract relevant data (e.g., "Genic", "Intergenic").

    Required Keys in `dataset`:
    - "{category} bases in SVs"
    - "{category} bases (genome)"
    - "Total bases in SVs"
    - "Total bases (genome)"
    """

    # Calculate proportions
    proportions = [dataset[1]/dataset[4],
                   dataset[0]/dataset[3],
                   dataset[2]/dataset[5]]

    non_proportions = [1 - p for p in proportions]

    # Plot
    plot_expected_vs_observed(ax, title, proportions, non_proportions)
    return proportions

def ghost_axis(ax, title, pad=20):
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.set_facecolor('none')
    ax.set_title(title, size=20, pad=pad,fontweight='bold')

def mod_hatches(ax):
    for i, bar in enumerate(ax.patches):
        if i == 0: # 0070bd
            bar.set_facecolor("#0070bd")  # Set color
            bar.set_hatch('')
        if i == 2: 
            bar.set_facecolor("#0070bd")  # Set color
            bar.set_hatch('///')
        if i == 1:
            bar.set_facecolor("#D2520D")  # Set color
            bar.set_hatch('')
        if i == 3:
            bar.set_facecolor("#D2520D")  # Set another color
            bar.set_hatch('///')       # Set another ha

def plot_scale(ax1_leg, ax1_leg1):
    color1 =  ["white", "#D2520D","#C94600"] 
    color2 = ["white", "#004FBD","#003887"]

    cmap = LinearSegmentedColormap.from_list("colormap", color1)
    cmap2 = LinearSegmentedColormap.from_list("colormap", color2)

    # Display the colormap
    gradient = np.linspace(0, 1, 6).reshape(1, -1)
    ax1_leg.imshow(gradient,  cmap=cmap2)
    ax1_leg.set_xticks([0, 5], ["low", "high"])
    ax1_leg.set_yticks([])
    ax1_leg.set_xlabel("$C. missioneira$", size=16)

    ax1_leg1.imshow(gradient,  cmap=cmap)
    ax1_leg1.set_xticks([0, 5], ["low", "high"])
    ax1_leg1.set_yticks([]) 
    ax1_leg1.set_xlabel("$C. tuca$", size=16)


# Define file paths and sample names
file_paths = [ variants, variants2]
samples = ["$C. tuca$", "$C. missioneira$"]

# Call the function
variants_plot = process_variants(file_paths, samples)

fig = plt.figure(figsize=(18, 20))
gs = GridSpec(4, 5, figure=fig, hspace=0.5)

variants_plot.groupby(["sample","type"]).describe()["length"]["mean"]
variants_plot.groupby(["sample","type"]).describe()["start"]["count"]


# subplots arrangement
ax1 = fig.add_subplot(gs[0, 0:4])
ax1_2 = fig.add_subplot(gs[1, 0:4])
ax1_leg = fig.add_subplot(gs[0:2,4])
ax1_divider = make_axes_locatable(ax1_leg)
ax1_leg1 = ax1_divider.append_axes("top", size="100%", pad="300%")

# ax0 = fig.add_subplot(gs[1, 4])
ax2 = fig.add_subplot(gs[3, 0:2])
ax3 = fig.add_subplot(gs[3, 2:4])
labels1 = fig.add_subplot(gs[3, 4]) 
ghost1 = fig.add_subplot(gs[3, 0:4])

ax4 = fig.add_subplot(gs[2, 0])
ax5 = fig.add_subplot(gs[2, 1])
ax6 = fig.add_subplot(gs[2, 2])
ax7 = fig.add_subplot(gs[2, 3])
labels = fig.add_subplot(gs[2, 4])
ghost2  = fig.add_subplot(gs[2, 1:4])


## start plotting

#add karyotypes
heatmap1 = np.asarray(Image.open(path + 'chromosome.png'))
heatmap = np.asarray(Image.open(path + 'chromosome2.png'))

ax1.imshow(heatmap1, aspect=0.30)
# ax1.set_xlim(right=2050, left=400)
# ax1.set_ylim(bottom=2005, top=600)
mod_hatches(ax1)
ax1.set_title("Normalized SV density", pad=30, size=20,fontweight='bold')
ax1.set_yticks([])
ax1.set_xticks(range(69, 1900,157), list(range(1,13)))

ax1_2.imshow(heatmap, aspect=0.30)
# ax1_2.set_xlim(right=2050, left=400)
# ax1_2.set_ylim(bottom=2005, top=600)
mod_hatches(ax1_2)
ax1_2.set_yticks([])
ax1_2.set_xticks(range(69, 1900,157), list(range(13,25)))

plot_scale(ax1_leg, ax1_leg1)

# plot sv count

fig1_df = variants_plot.groupby(["sample", "genotype"]).size().unstack()
fig1_df[["1/1","0/1"]].plot.bar(width=0.5, ax=ax4)
patterns = ['/', '/']  # Define patterns for each group
mod_hatches(ax4)

###### 2nd line plots per species ######
###### sv proportions 
######## intergenic

intergenic = [Ctuca_results["Intergenic bases in SVs"], 
            C_missioneira_results["Intergenic bases in SVs"],
            Ctuca_results["Intergenic bases (genome)"],
            Ctuca_results["Total SV bases"],
            C_missioneira_results["Total SV bases"],
            Ctuca_results["Total genome bases"]]


p = plot_proportions(ax5, "Intergenic | Genic", intergenic)
ax5.legend().set_visible(False)

######## genic
# genic = [Ctuca_results["Genic bases in SVs"], 
#          C_missioneira_results["Genic bases in SVs"],
#          Ctuca_results["Genic bases (genome)"],
#          Ctuca_results["Total SV bases"],
#          C_missioneira_results["Total SV bases"],
#          Ctuca_results["Total genome bases"]]

# plot_proportions(ax5, "Genic", genic)

######## intronic
intronic = [Ctuca_results2["Intronic bases in SVs"], 
            C_missioneira_results2["Intronic bases in SVs"],
            Ctuca_results2["Intronic bases (genome)"],
            Ctuca_results["Genic bases in SVs"],
            C_missioneira_results["Genic bases in SVs"],
            Ctuca_results["Genic bases (genome)"]]

plot_proportions(ax6, "Intronic | Exonic", intronic)
ax6.set_yticks([])

# ######## exonic

# exonic = [Ctuca_results2["Exonic bases in SVs"],
#           C_missioneira_results2["Exonic bases in SVs"],
#             Ctuca_results2["Exonic bases (genome)"],
#             Ctuca_results["Genic bases in SVs"],
#             C_missioneira_results["Genic bases in SVs"],
#              Ctuca_results["Genic bases (genome)"]]

# plot_proportions(ax7, "Exonic", exonic)

#### repeats

repeats = [Ctuca_results3["Repeat bases in SVs"],
           C_missioneira_results3["Repeat bases in SVs"],
           repeat_bases_in_genome, 
           Ctuca_results["Total SV bases"],
           C_missioneira_results["Total SV bases"],
           total_genome_bases]
plot_proportions(ax7, "Repeats | Non-repeats", repeats)
ax7.legend().remove()
### plot SV length distribution
categorical = variants_plot.groupby(["sample", "length_category", "type"]).size().unstack().reset_index()

plot_stacked_prop(categorical, "$C. missioneira$", ax2)
plot_stacked_prop(categorical, "$C. tuca$", ax3)

ghost_axis(ghost1, "SV lenght distribution", 30)
ghost_axis(ghost2, "SV proportions", 30)

########### legends and labels
# remove al innecesary lines
format_axes([ax1,ax1_2,ax2,ax3,ax4,ax5,ax6,ax7,labels, labels1])
ax4.set_yticks([0,50000, 100000, 150000], ["0","50k","100k", "150k"])
# ax4.yaxis.tick_right()
patch1 = mpatches.Patch(color="#ff5900", label='$C. missioneira$')
patch3 = mpatches.Patch(hatch='/', label='0/1', edgecolor='black', facecolor='darkgray')
patch2 = mpatches.Patch(color="#D2520D", label='$C. tuca$')
patch4 = mpatches.Patch(facecolor='darkgray', label='1/1',edgecolor='black')
ax4.legend(handles=[patch3,patch4], ncols=1,frameon=False, fancybox=False, fontsize=16, bbox_to_anchor=(0.3, 0.9))
ax4.set_title("SV count", size=20,fontweight='bold', pad=30)

ax2.legend().remove()
ax2.set_title("$C. missioneira$", size=16)
ax2.yaxis.set_visible(False)

ax3.set_title("$C. tuca$",size=16)
ax3.yaxis.tick_right()
ax3.set_ylabel("")
ax3.spines['right'].set_position(('axes', 0.95))
ax3.set_ymargin(0)
leg = ax3.legend(ncols=1,frameon=False, fancybox=False, fontsize=16)
labels1.legend(leg.legendHandles, ["Duplication","Indel","Inversion"], ncols=1,frameon=False, fancybox=False, fontsize=16, loc = "upper right") 
ax3.legend().remove()
patch1 = mpatches.Patch(color="#0070bd", label='$C. missioneira$')
patch2 = mpatches.Patch(color="#D2520D", label='$C. tuca$')
patch3 = mpatches.Patch(color='darkgray', label='Background')
labels.legend(handles=[patch1, patch2, patch3], ncols=1,frameon=False, fancybox=False, fontsize=16, loc = "upper right")


ax7.set_yticks([0, 0.5, 1])
ax7.yaxis.tick_right()

labels.yaxis.set_visible(False)
labels.xaxis.set_visible(False)
labels1.yaxis.set_visible(False)
labels1.xaxis.set_visible(False)
ax5.legend().remove()
ax6.legend().remove()
ax7.legend().remove()
labels1.set_facecolor('none')


# plt.savefig(path + "SV_genome_report.svg", dpi=300)

# plt.savefig(path + "SV_genome_report.png", dpi=300)
