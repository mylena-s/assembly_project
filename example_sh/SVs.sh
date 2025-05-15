# # C. tuca
./nextflow run StructuralVariation.nf --genome CBFH00195_synteny.fasta --reads data/CBFH00195_PAQ76882.fastq_trimmed.fq.gz  --publishDir results  --covmax 45 --covmin 15 --sniffles
# # C. missioneira
./nextflow run StructuralVariation.nf --genome CBFH00195_synteny.fasta --reads data/CBFH00251_PAQ76522.fastq_trimmed.fq.gz --publishDir results --covmax 30 --covmin 10 --sniffles

#run again after directory error for sniffles
# C. tuca
#./nextflow run StructuralVariation.nf --genome CBFH00195_synteny.fasta --mapped results/alignment/CBFH00195_synteny_CBFH00251_PAQ76522_trimmed.fq_sorted.bam  --publishDir results --debreak --covmax 40 --covmin 20 --resume
# C. missioneira
#./nextflow run StructuralVariation.nf --genome CBFH00195_synteny.fasta --mapped results/alignment/CBFH00195_synteny_CBFH00195_PAQ76882_trimmed.fq_sorted.bam --publishDir results --debreak --covmax 30 --covmin 10


# run all changing to winnowmap aligner
#./nextflow run StructuralVariation.nf --genome CBFH00195_synteny.fasta --reads data/CBFH00195_PAQ76882.fastq_trimmed.fq.gz  --publishDir results --covmax 40 --covmin 20
#./nextflow run StructuralVariation.nf --genome CBFH00195_synteny.fasta --reads data/CBFH00251_PAQ76522.fastq_trimmed.fq.gz --publishDir results --covmax 30 --covmin 10


# run from alignment
#./nextflow run StructuralVariation.nf --genome CBFH00195_synteny.fasta --mapped results/alignment/CBFH00195_synteny_CBFH00195_PAQ76882_winnowmap_sorted.bam --publishDir results --covmax 40 --covmin 20 --sniffles
#./nextflow run StructuralVariation.nf --genome CBFH00195_synteny.fasta --mapped results/alignment/CBFH00195_synteny_CBFH00251_PAQ76522_winnowmap_sorted.bam --publishDir results --covmax 30 --covmin 10 --debreak --sniffles
