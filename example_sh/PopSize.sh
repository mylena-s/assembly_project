# usando scripts de msmc
#./nextflow run PopSize.nf --genome CBFH00195_renamed.fasta --bam results/alignment/CBFH00195_renamed_CBFH00195_trimmed.fq_sorted.bam --satellites results/centromere/CBFH00195_satellite_renamed2.bed --publishDir results

# usando deepvariant vcf filtrado para calidad y heterosigosis
./nextflow run PopSize.nf --genome CBFH00195_renamed.fasta --bam results/alignment/CBFH00195.fastq_trimmed.sorted.q30.bam --satellites results/centromere/CBFH00195_satellite_renamed.bed --vcf /media/jmf/DATOS2/mylena/genome_report/results/population/heterozygosis/CBFH00195_all_joined_filtered.het.vcf.gz --publishDir results --meandep 34  --contigs CBFH00195_synteny.contigs