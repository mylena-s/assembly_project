
process CENTROMERE{
    label 'medium_resources'
    publishDir "${params.publishDir}/annotation/repeat/", mode: 'copy', pattern: '03_nessie*'
    errorStrategy 'ignore'
    label 'minimap'
    input:
    path genome
    each mode

    output:
    path "03_nessie_${genome.baseName}"
    path "03_nessie_${genome.baseName}/*_complexity.tsv", optional: true, emit: complexity
    path "03_nessie_${genome.baseName}/*_entropy.tsv", optional: true, emit: entropy

    script:
    """
    mkdir 03_nessie_${genome.baseName}
    nessie -I $genome -O $genome$mode-l 10000 -s 10000
    for file in *tsv; do to_wig.py -i \$file -o \$file'.wig'; done
    for file in *wig; do tail -n +2 \$file | grep -v 'track type' - > \$file'.bw'; done
    samtools faidx $genome
    awk '{print \$1, \$2}' $genome'.fai' > $genome'.chrsizes'
    for file in *.bw; do wigToBigWig \$file *.chrsizes \$file'.bigwig'; done
    for file in *bigwig; do bigWigToBedGraph \$file \$file'.bedgraph'; done
    cp *tsv 03_nessie_${genome.baseName}
    cp *bigwig 03_nessie_${genome.baseName}
    cp *bedgraph 03_nessie_${genome.baseName}
    echo "statistics in bedgraph correspond to the midpoint of the genomic windows" >> .command.log
    cp .command.sh 03_nessie_${genome.baseName}/.${genome.baseName}.command.sh
    cp .command.log 03_nessie_${genome.baseName}/.${genome.baseName}.command.log
    """
}
