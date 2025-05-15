./nextflow run Ribosomal.nf --ref /media/jmf/DATOS2/mylena/genome_report/results/rna/exconsensus.fa \
    --hifiasm_folder /media/jmf/DATOS2/mylena/genome_report/results/rna/hifiasm \
    --hifiasm_prefix CBFH00195.fastq_trimmed.fq.asm\
    --tangle /media/jmf/DATOS2/mylena/genome_report/results/rna/tangle1.txt \
    --ont /media/jmf/DATOS2/mylena/genome_report/results/rna/CBFH00195_PAQ76882.fastq_trimmed.fasta \
    --hifi /media/jmf/DATOS2/mylena/genome_report/results/rna/CBFH00195.fastq_trimmed.fasta \
    --publishDir results --fasta true
