# gene annotation after repeat masking with braker2
./nextflow Annotation.nf -c nextflow_ferropc.config --publishDir CBFH00195/results/ --genome CBFH00195/results/annotation/repeat/04_masking/CBFH00195.fasta.masked --masked true 
# gene annotation with galba
./nextflow Annotation.nf -c nextflow_ferropc.config --publishDir CBFH00195/results/ --genome CBFH00195/results/annotation/repeat/04_masking/CBFH00195.fasta.masked --masked true --proteins lib/Cichliformes.fasta.gz --mode galba

# functional gene annotation after galba
./nextflow Annotation.nf -c nextflow_ferropc.config --publishDir CBFH00195/results/ --genome CBFH00195/results/annotation/repeat/04_masking/CBFH00195.fasta.masked --masked true --protannot CBFH00195/results/annotation/gene/structural/GALBA/galba.aa

# functional annotation after TOGA and busco run
./nextflow Annotation.nf -c nextflow_ferropc.config --publishDir CBFH00195/results/ --genome CBFH00195/results/annotation/repeat/04_masking/CBFH00195.fasta.masked --masked true --protannot CBFH00195/results/annotation/gene/toga/*fasta
