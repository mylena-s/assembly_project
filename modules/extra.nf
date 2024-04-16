process ALIGN_CACTUS { 
    label 'cactus'
    label 'medium_resources'
    
    input:
    path genomes

    output:
    path "mitogenomes.hal"
    path "mitogenomes.maf.gz"
    
    script:
    def files = genomes.join(' ')
    """
    ls $files > genomes.txt
    cactus ./js ./genomes.txt ./mitogenomes.hal --maxMemory $task.memory
    cactus-hal2maf ./js mitogenomes.hal mitogenomes.maf.gz --refGenome CBFH00195 --chunkSize 1000000"""
}


process ALIGN_GENOMES {
    label 'medium_resources'
    label 'minimap'
    
    input:
    path genomes

    output:
    path alignment

    script:
    def files = genomes.join(' ')
    """
    cat $files > genomes.fasta
    minimap2 -ax asm5 CBFH00195.final_mitogenome.fasta genomes.fasta > multiple_mito.sam
    """
}