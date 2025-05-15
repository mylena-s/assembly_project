# complete run
# ./nextflow run HaplotypePhasing.nf -c nextflow.config --reference CBFH00195_renamed.fasta --resume --hifi data/CBFH00195.fastq_trimmed.fq.gz --ont data/CBFH00195_PAQ76882.fastq_trimmed.fq.gz --publishDir results 

# run after variant calling
./nextflow run Phylogeny.nf -c nextflow.config --genome CBFH00195_synteny.fasta --publishDir results \
    --vcf results/population/glnexus/ont_allmerged.vcf.gz \
    --bams "results/population/glnexus/*bam" --twiss no --phasing hifi