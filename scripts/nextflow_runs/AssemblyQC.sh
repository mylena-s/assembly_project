# RUN ALL STEPS
./nextflow run AssemblyQC.nf --publishDir CBFH00195/results --reads CBFH00195/results/pre_assembly/trimming/CBFH00195_PAQ76882.fastq_trimmed.fq.gz -c nextflow.config --genome "/home/fhenning/assembly_project/CBFH00195/results/assembly/QC/assemblies/*.fasta" --all -resume
# RUN ONLY CORRECTNESS ESTIMATES
./nextflow run AssemblyQC.nf --publishDir CBFH00195/results --reads CBFH00195/results/pre_assembly/trimming/CBFH00195_PAQ76882.fastq_trimmed.fq.gz -c nextflow.config --genome "/home/fhenning/assembly_project/CBFH00195/results/assembly/QC/assemblies/*.fasta" --correctness
# TUN ONLY CONTAMINATION ESTIMATES
./nextflow run AssemblyQC.nf --publishDir CBFH00195/results --reads CBFH00195/results/pre_assembly/trimming/CBFH00195_PAQ76882.fastq_trimmed.fq.gz -c nextflow.config --genome "/home/fhenning/assembly_project/CBFH00195/results/assembly/QC/assemblies/*.fasta" --cleaness
