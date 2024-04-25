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

process SAMPLE_ONT {
    label 'low_resources'
    label 'nanoplot'

    input:
    path reads
    val (lreads)
    
    output:
    path "*subset.fastq.gz", emit:reads
    script:
    """gunzip -c $reads | NanoFilt -l $lreads | pigz > ${reads.baseName}_${lreads}_subset.fastq.gz"""

}


process SCAFFOLDING {
    label 'goldrush'
    label 'medium_resources'
    publishDir "${params.publishDir}/assembly/ntlink_scaffolding/", mode: 'copy'  
    errorStrategy 'ignore'

    input:
    each path(genome) 
    path reads,  name: "reads.fq.gz"

    output:
    path("ntlink_${genome}"), emit: scaffolded_genome
    path "ntlink_${genome}.agp", emit: assemblygraph
    path ".*command*"

    script:
    """
    seqkit seq -m 1000 $genome > genome_subset.fasta
    ntLink_rounds run_rounds_gaps target=genome_subset.fasta reads=reads.fq.gz k=32 w=250 t=$task.cpus rounds=4 clean=True overlap=True a=$params.nreads paf=True
    mv genome_subset.fasta.k32.w250.z1000.ntLink.ntLink.gap_fill.fa.k32.w250.z1000.ntLink.scaffolds.gap_fill.fa ntlink_${genome}
    mv genome_subset.fasta.k32.w250.z1000.ntLink.ntLink.gap_fill.fa.k32.w250.z1000.ntLink.scaffolds.gap_fill.fa.agp ntlink_${genome}.agp
    mv *paf ntlink_${genome}.paf
    mv .command.sh .${genome}.command.sh
    mv .command.log .${genome}.command.log
    """
}

params.nreads = 2
