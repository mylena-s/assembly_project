# whole dataset, synteny figure
./nextflow run Synteny.nf --reference CBFH00195_synteny.fasta --query Amphilopus_synteny.fasta --query2 Oniloticus_synteny.fasta --publishDir results --resume --genomes_file genomes.txt
# all chr except chr10, for matrix estimation
./nextflow run Synteny.nf --reference CBFH00195_synteny_no10.fasta --query2 Amphilopus_synteny_no10.fasta --query Oniloticus_synteny_no10.fasta --publishDir results --plot no --genomes_file genomes.txt
./nextflow run Synteny.nf --reference CBFH00195_synteny_no10.fasta --query Amphilopus_synteny_no10.fasta --query2 Oniloticus_synteny_no10.fasta --publishDir results --plot no --genomes_file genomes.txt
