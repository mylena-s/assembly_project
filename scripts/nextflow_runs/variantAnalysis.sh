# # missioneira run
# ./nextflow run Variant_analysis.nf --genome CBFH00195_synteny.fasta --hom --call --variant_analysis --bam results/alignment/CBFH00195_synteny_CBFH00251_PAQ76522_trimmed.fq_sorted.bam --publishDir results --maxDepth 30 --minDepth 7 --contigs CBFH00195_synteny.list

## missionaira from calling
# ./nextflow run Variant_analysis.nf --genome CBFH00195_synteny.fasta --hom --variant_analysis --popsize --resume --bam results/alignment/CBFH00195_synteny_CBFH00251_PAQ76522_trimmed.fq_sorted.bam --vcf results/population/variant_calling/CBFH00251.all_joined.vcf.gz --publishDir results --maxDepth 30 --minDepth 7 --contigs CBFH00195_synteny.list --satellites /media/jmf/DATOS2/mylena/genome_report/results/centromere/CBFH00195_synteny.srf-aln.bed


# # tuca ont, call 
#./nextflow run Variant_analysis.nf --genome CBFH00195_synteny.fasta --hom --call --bam results/alignment/CBFH00195_synteny_CBFH00195_PAQ76882_trimmed.fq_sorted.bam --publishDir results --maxDepth 40 --minDepth 15 --contigs CBFH00195_synteny.list

# tuca ont, from variant calling to pop estimates
#./nextflow run Variant_analysis.nf --genome CBFH00195_synteny.fasta --hom --popsize --variant_analysis --bam results/alignment/CBFH00195_synteny_CBFH00195_PAQ76882_trimmed.fq_sorted.bam --vcf /media/jmf/DATOS2/mylena/genome_report/results/population/variant_calling/CBFH00195_all_joined.vcf.gz --publishDir results --maxDepth 40 --minDepth 15 --contigs CBFH00195_synteny.list --satellites /media/jmf/DATOS2/mylena/genome_report/results/centromere/CBFH00195_synteny.srf-aln.bed


# tuca hifi 
# ./nextflow run Variant_analysis.nf --genome CBFH00195_synteny.fasta --reads data/CBFH00195.fastq_trimmed.fq.gz --popsize --variant_analysis --hom --call --publishDir results --maxDepth 50 --minDepth 10 --contigs CBFH00195_synteny.list --satellites /media/jmf/DATOS2/mylena/genome_report/results/centromere/CBFH00195_synteny.srf-aln.bed

# tuca hifi  from variant calling to popsize
#./nextflow run Variant_analysis.nf --genome CBFH00195_synteny.fasta --popsize --variant_analysis --hom --variant_analysis --bam results/alignment/CBFH00195_synteny_CBFH00195_trimmed.fq_sorted.bam --vcf /media/jmf/DATOS2/mylena/genome_report/results/population/variant_calling/CBFH00195_synteny_CBFH00195_trimmed_all_joined.vcf.gz --publishDir results --maxDepth 50 --minDepth 10 --contigs CBFH00195_synteny.list --satellites /media/jmf/DATOS2/mylena/genome_report/results/centromere/CBFH00195_synteny.srf-aln.bed

