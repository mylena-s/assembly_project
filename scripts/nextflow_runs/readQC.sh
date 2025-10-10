export PATH="/home/fhenning/Software":$PATH

./nextflow run readQC.nf --publishDir CBFH00195/results --reads CBFH00195/data/CBFH00195_PAQ76882.fastq.gz -c nextflow.config
