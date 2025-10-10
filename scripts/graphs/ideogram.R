library("RIdeogram")
library(dplyr)
library(stringr)
library(preprocessCore)

setwd('/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/')

# normalize sv density

SVSbed195 = read.table("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/processed/CBFH00195_coverage.bed", sep="\t", header=FALSE, colClasses = c("character", "numeric",
    "numeric","numeric","numeric","numeric", "numeric"))

SVSbed251 = read.table("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/processed/CBFH00251_coverage.bed", sep="\t", header=FALSE, colClasses = c("character", "numeric",
    "numeric","numeric","numeric","numeric","numeric"))
INM = read.table("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/tracks/toga.imm.genes.cov.bed",sep="\t", header=FALSE, colClasses = c("character", "numeric", "numeric","numeric","numeric","numeric","numeric"))

colnames(INM) = c("Chr", "Start", "End","X","XX","XXX","Value")
INM = INM %>% mutate(across('Chr', str_replace, 'LG', ''))
INM$Chr = as.numeric(INM$Chr)
INM = INM %>% arrange(Chr)


SVSbed195 = SVSbed195 %>% select(V1,V2,V3,V4)
SVSbed251 = SVSbed251 %>% select(V1,V2,V3,V4)
head(SVSbed251)
head(SVSbed195)

df <- merge(SVSbed195, SVSbed251, by = c("V1", "V2", "V3"), suffixes = c("_195", "_251"))
df = df[order(df$V1,df$V2),]
str(df)
norm = normalize.quantiles(as.matrix(df[,c("V4_195","V4_251")]))

merged_df = cbind(df, norm)

colnames(merged_df) = c("chr", "start", "end", "195", "251", "195_norm", "251_norm")


karyo = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/karyotype/CBFH00195.chrominfo"
chrms <- read.table(karyo, sep = '\t', header = F, stringsAsFactors = F)
centromere <- read.table("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/karyotype/centromere.bed", sep = '\t', header = T, stringsAsFactors = F)


colnames(chrms) = c('Chr', 'Start', 'End')
chrms = chrms %>% mutate(across('Chr', str_replace, 'LG', ''))
chrms =  merge(chrms, centromere, by = "Chr", all=T)

svdensity <- merged_df
str(svdensity)

svdensity = svdensity %>% mutate(across('chr', str_replace, 'LG', ''))
svdensity$chr = as.numeric(svdensity$chr)


Ctuca = svdensity %>% select(chr, start, end, '195_norm')
Cmissioneira = svdensity %>% select(chr, start, end, '251_norm')

colnames(Ctuca) = c('Chr', 'Start', 'End', 'Value')
colnames(Cmissioneira) = c('Chr', 'Start', 'End', 'Value')

chrms$Chr = as.numeric(chrms$Chr)
chrms = chrms %>% arrange(Chr)


ideogram(karyotype = chrms[13:24,], overlaid = Ctuca, label = Cmissioneira, label_type = "heatmap",  colorset1 = c("white", "#D2520D","#C94600","black"), colorset2 = c("white", "#004FBD","#003887","black"), Ly=-2, Lx=-2)
# convertSVG("chromosome.svg", device="png", dpi = 300)

#change name and run 

ideogram(karyotype = chrms[1:12,], overlaid = Ctuca, label = Cmissioneira, label_type = "heatmap", colorset1 = c("white", "#D2520D","#C94600","black"), colorset2 = c("white", "#004FBD","#003887","black"), Ly = 1)
# convertSVG("chromosome.svg", device="png", dpi = 300)


### with markers
#ARE INMUNOGLOBULINGS THE CAUSE OF MISSASSEMBLY INTERPRETED AS SVS? NO
INM = INM %>% select(Chr, Start, End, Value)
ideogram(karyotype = chrms[1:12,], overlaid = Ctuca, label = INM, label_type = "heatmap", colorset1 = c("white", "#D2520D","#C94600","black"), colorset2 = c("white", "#004FBD","#003887","black"), Ly = 1)
ideogram(karyotype = chrms[13:24,], overlaid = Ctuca, label = INM, label_type = "heatmap",  colorset1 = c("white", "#D2520D","#C94600","black"), colorset2 = c("white", "#004FBD","#003887","black"), Ly=-2, Lx=-2)


#### TRANS SPECIES POLYMORPHISM IN INMUNOGLOBULIN GENES


TSP = read.table("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/processed/suppvec11_coverage.bed",sep="\t", header=FALSE, colClasses = c("character", "numeric", "numeric","numeric","numeric","numeric","numeric"))
colnames(TSP) = c("Chr", "Start", "End","X","XX","XXX","Value")
TSP = TSP %>% mutate(across('Chr', str_replace, 'LG', ''))
TSP$Chr = as.numeric(TSP$Chr)
TSP = TSP %>% arrange(Chr)

ideogram(karyotype = chrms[1:12,], overlaid = INM, label = TSP, label_type = "heatmap", colorset1 = c("white", "#D2520D","#C94600","black"), colorset2 = c("white", "#004FBD","#003887","black"), Ly = 1)
ideogram(karyotype = chrms[13:24,], overlaid = INM, label = TSP, label_type = "heatmap",  colorset1 = c("white", "#D2520D","#C94600","black"), colorset2 = c("white", "#004FBD","#003887","black"), Ly=-2, Lx=-2)

### SVs and repeats
REP = read.table("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/tracks/repeats.final.coverage.bed",sep="\t", header=FALSE, colClasses = c("character", "numeric", "numeric","numeric","numeric","numeric","numeric"))
colnames(REP) = c("Chr", "Start", "End","X","XX","XXX","Value")
REP = REP %>% mutate(across('Chr', str_replace, 'LG', ''))
REP$Chr = as.numeric(REP$Chr)
REP = REP %>% arrange(Chr)
ideogram(karyotype = chrms[1:12,], overlaid = Ctuca, label = REP, label_type = "heatmap", colorset1 = c("white", "#D2520D","#C94600","black"), colorset2 = c("white", "#2CA25F", "#006D2C", "black"), Ly = 1)
ideogram(karyotype = chrms[13:24,], overlaid = Ctuca, label = REP, label_type = "heatmap",  colorset1 = c("white", "#D2520D","#C94600","black"), colorset2 = c("white", "#2CA25F", "#006D2C", "black"), Ly=-2, Lx=-2)

### SVs and sats
REP = read.table("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/repeats/satellite.bed.cov.bed",sep="\t", header=FALSE, colClasses = c("character", "numeric", "numeric","numeric","numeric","numeric","numeric"))
colnames(REP) = c("Chr", "Start", "End","X","XX","XXX","Value")
REP = REP %>% mutate(across('Chr', str_replace, 'LG', ''))
REP$Chr = as.numeric(REP$Chr)
REP = REP %>% arrange(Chr)
ideogram(karyotype = chrms[1:12,], overlaid = Ctuca, label = REP, label_type = "heatmap", colorset1 = c("white", "#D2520D","#C94600","black"), colorset2 = c("white", "#2CA25F", "#006D2C", "black"), Ly = 1)
ideogram(karyotype = chrms[13:24,], overlaid = Ctuca, label = REP, label_type = "heatmap",  colorset1 = c("white", "#D2520D","#C94600","black"), colorset2 = c("white", "#2CA25F", "#006D2C", "black"), Ly=-2, Lx=-2)

### SVs and LINES
REP = read.table("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/repeats/repeats.all.coverage.bed",sep="\t", header=FALSE, colClasses = c("character", "numeric", "numeric","numeric","numeric","numeric","numeric"))
colnames(REP) = c("Chr", "Start", "End","X","XX","XXX","Value")
REP = REP %>% mutate(across('Chr', str_replace, 'LG', ''))
REP$Chr = as.numeric(REP$Chr)
REP = REP %>% arrange(Chr)
ideogram(karyotype = chrms[1:12,], overlaid = Ctuca, label = REP, label_type = "heatmap", colorset1 = c("white", "#D2520D","#C94600","black"), colorset2 = c("white", "#2CA25F", "#006D2C", "black"), Ly = 1)
ideogram(karyotype = chrms[13:24,], overlaid = Ctuca, label = REP, label_type = "heatmap",  colorset1 = c("white", "#D2520D","#C94600","black"), colorset2 = c("white", "#2CA25F", "#006D2C", "black"), Ly=-2, Lx=-2)

############ CHR 3 + 23 
REP = read.table("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/repeats/DNA.bed.cov.bed",sep="\t", header=FALSE, colClasses = c("character", "numeric", "numeric","numeric","numeric","numeric","numeric"))
colnames(REP) = c("Chr", "Start", "End","X","XX","XXX","Value")
REP = REP %>% mutate(across('Chr', str_replace, 'LG', ''))
REP$Chr = as.numeric(REP$Chr)
REP = REP %>% arrange(Chr)
ideogram(karyotype = chrms[c(3,23),], overlaid = Ctuca, label = REP, label_type = "heatmap", colorset1 =c("white","black"), colorset2 =c("white", "#D2520D","#C94600","black"), Ly = 1, width = 50)





ideogram(karyotype = chrms[c(3,23),], overlaid = Ctuca, label = INM, label_type = "heatmap", colorset1 =c("white","black"), colorset2 =c("white", "#D2520D","#C94600","black"), Ly = 1, width = 50)
ideogram(karyotype = chrms[c(3,23),], overlaid = Ctuca, label = TSP, label_type = "heatmap", colorset1 =c("white","black"), colorset2 =c("white", "#D2520D","#C94600","black"), Ly = 1, width = 50)

########### SNPS


setwd('/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/')

# normalize sv density

SVSbed195 = read.table("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/heterozygosis/CBFH00195_all_joined_filtered.het.bed", sep="\t", header=FALSE, colClasses = c("character", "numeric",
    "numeric","numeric","numeric","numeric", "numeric"))

SVSbed251 = read.table("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/heterozygosis/CBFH00251_filtered.het.bed", sep="\t", header=FALSE, colClasses = c("character", "numeric",
    "numeric","numeric","numeric","numeric","numeric"))

SVSbed195 = SVSbed195 %>% select(V1,V2,V3,V4)
SVSbed251 = SVSbed251 %>% select(V1,V2,V3,V4)
head(SVSbed251)
head(SVSbed195)

df <- merge(SVSbed195, SVSbed251, by = c("V1", "V2", "V3"), suffixes = c("_195", "_251"))
df = df[order(df$V1,df$V2),]
str(df)
norm = normalize.quantiles(as.matrix(df[,c("V4_195","V4_251")]))

merged_df = cbind(df, norm)

colnames(merged_df) = c("chr", "start", "end", "195", "251", "195_norm", "251_norm")


karyo = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/karyotype/CBFH00195.chrominfo"
chrms <- read.table(karyo, sep = '\t', header = F, stringsAsFactors = F)
centromere <- read.table("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/karyotype/centromere.bed", sep = '\t', header = T, stringsAsFactors = F)


colnames(chrms) = c('Chr', 'Start', 'End')
chrms = chrms %>% mutate(across('Chr', str_replace, 'LG', ''))
chrms =  merge(chrms, centromere, by = "Chr", all=T)

svdensity <- merged_df
str(svdensity)

svdensity = svdensity %>% mutate(across('chr', str_replace, 'LG', ''))
svdensity$chr = as.numeric(svdensity$chr)

Ctuca = svdensity %>% select(chr, start, end, '195_norm')
Cmissioneira = svdensity %>% select(chr, start, end, '251_norm')

colnames(Ctuca) = c('Chr', 'Start', 'End', 'Value')
colnames(Cmissioneira) = c('Chr', 'Start', 'End', 'Value')

chrms$Chr = as.numeric(chrms$Chr)
chrms = chrms %>% arrange(Chr)


ideogram(karyotype = chrms[13:24,], overlaid = Ctuca, label = Cmissioneira, label_type = "heatmap",  colorset1 = c("white", "#D2520D","#C94600","black"), colorset2 = c("white", "#004FBD","#003887","black"), Ly=-2, Lx=-2)
convertSVG("chromosome.svg", device="png", dpi = 300)

#change name and run 

ideogram(karyotype = chrms[1:12,], overlaid = Ctuca, label = Cmissioneira, label_type = "heatmap", colorset1 = c("white", "#D2520D","#C94600","black"), colorset2 = c("white", "#004FBD","#003887","black"), Ly = 1)
convertSVG("chromosome.svg", device="png", dpi = 300)



########## snps & Svs

########## tuca

setwd('/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/')

# normalize sv density

SVSbed195 = read.table("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/processed/CBFH00195_coverage.bed", sep="\t", header=FALSE, colClasses = c("character", "numeric",
    "numeric","numeric","numeric","numeric", "numeric"))

SNPsbed195 = read.table("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/heterozygosis/CBFH00195_all_joined_filtered.het.bed", sep="\t", header=FALSE, colClasses = c("character", "numeric",
    "numeric","numeric","numeric","numeric","numeric"))

SVSbed195 = SVSbed195 %>% select(V1,V2,V3,V4)
SNPsbed195 = SNPsbed195 %>% select(V1,V2,V3,V4)
head(SNPsbed195)
head(SVSbed195)

df <- merge(SVSbed195, SNPsbed195, by = c("V1", "V2", "V3"), suffixes = c("_SV", "_SNP"))
df = df[order(df$V1,df$V2),]
str(df)
norm = normalize.quantiles(as.matrix(df[,c("V4_SV","V4_SNP")]))

merged_df = cbind(df, norm)

colnames(merged_df) = c("chr", "start", "end", "195", "251", "SV_norm", "SNP_norm")


karyo = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/karyotype/CBFH00195.chrominfo"
chrms <- read.table(karyo, sep = '\t', header = F, stringsAsFactors = F)
centromere <- read.table("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/karyotype/centromere.bed", sep = '\t', header = T, stringsAsFactors = F)


colnames(chrms) = c('Chr', 'Start', 'End')
chrms = chrms %>% mutate(across('Chr', str_replace, 'LG', ''))
chrms =  merge(chrms, centromere, by = "Chr", all=T)

svdensity <- merged_df
str(svdensity)

svdensity = svdensity %>% mutate(across('chr', str_replace, 'LG', ''))
svdensity$chr = as.numeric(svdensity$chr)

SV = svdensity %>% select(chr, start, end, 'SV_norm')
SNP = svdensity %>% select(chr, start, end, 'SNP_norm')

colnames(SV) = c('Chr', 'Start', 'End', 'Value')
colnames(SNP) = c('Chr', 'Start', 'End', 'Value')

chrms$Chr = as.numeric(chrms$Chr)
chrms = chrms %>% arrange(Chr)


ideogram(karyotype = chrms[13:24,], overlaid = SV, label = SNP, label_type = "heatmap",  colorset1 = c("white", "#D2520D","#C94600","black"), colorset2 = c("white", "#004FBD","#003887","black"), Ly=-2, Lx=-2)
ideogram(karyotype = chrms[1:12,], overlaid = SV, label = SNP, label_type = "heatmap", colorset1 = c("white", "#D2520D","#C94600","black"), colorset2 = c("white", "#004FBD","#003887","black"), Ly = 1)



############################# missio
SVSbed251 = read.table("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/processed/CBFH00251_coverage.bed", sep="\t", header=FALSE, colClasses = c("character", "numeric",
    "numeric","numeric","numeric","numeric","numeric"))

SNPbed251 = read.table("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/heterozygosis/CBFH00251_filtered.het.bed", sep="\t", header=FALSE, colClasses = c("character", "numeric",
    "numeric","numeric","numeric","numeric","numeric"))

SVSbed251 = SVSbed251 %>% select(V1,V2,V3,V4)
SNPbed251 = SNPbed251 %>% select(V1,V2,V3,V4)


df <- merge(SVSbed251, SNPbed251, by = c("V1", "V2", "V3"), suffixes = c("_SV", "_SNP"))
df = df[order(df$V1,df$V2),]
str(df)
norm = normalize.quantiles(as.matrix(df[,c("V4_SV","V4_SNP")]))

merged_df = cbind(df, norm)

colnames(merged_df) = c("chr", "start", "end", "195", "251", "SV_norm", "SNP_norm")


karyo = "/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/karyotype/CBFH00195.chrominfo"
chrms <- read.table(karyo, sep = '\t', header = F, stringsAsFactors = F)
centromere <- read.table("/home/mylequis/Documents/Doctorado/Papers/Genome_report/00_DATA/Figures_ideograms/karyotype/centromere.bed", sep = '\t', header = T, stringsAsFactors = F)


colnames(chrms) = c('Chr', 'Start', 'End')
chrms = chrms %>% mutate(across('Chr', str_replace, 'LG', ''))
chrms =  merge(chrms, centromere, by = "Chr", all=T)

svdensity <- merged_df
str(svdensity)

svdensity = svdensity %>% mutate(across('chr', str_replace, 'LG', ''))
svdensity$chr = as.numeric(svdensity$chr)

SV = svdensity %>% select(chr, start, end, 'SV_norm')
SNP = svdensity %>% select(chr, start, end, 'SNP_norm')

colnames(SV) = c('Chr', 'Start', 'End', 'Value')
colnames(SNP) = c('Chr', 'Start', 'End', 'Value')

chrms$Chr = as.numeric(chrms$Chr)
chrms = chrms %>% arrange(Chr)


ideogram(karyotype = chrms[13:24,], overlaid = SV, label = SNP, label_type = "heatmap",  colorset1 = c("white", "#D2520D","#C94600","black"), colorset2 = c("white", "#004FBD","#003887","black"), Ly=-2, Lx=-2)
ideogram(karyotype = chrms[1:12,], overlaid = SV, label = SNP, label_type = "heatmap", colorset1 = c("white", "#D2520D","#C94600","black"), colorset2 = c("white", "#004FBD","#003887","black"), Ly = 1)
