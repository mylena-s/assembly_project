process MINIMAP2{
    label 'medium_resources'
    label 'syri'
    maxForks 2
    publishDir "${params.publishDir}/alignment", mode: 'copy'   


    input:
    tuple path(genome1), path(genome2)
    
    output:
    tuple path ('*.bam'), path("*bai"), emit: mappings
    path '.*command.*'
    script:
    def name = "${genome1.simpleName}_${genome2.simpleName}"

    """minimap2 -ax asm5 -f 0.02 -t 10 --eqx $genome1 $genome2 | samtools sort -O BAM - > ${name}.bam
    samtools index ${name}.bam
    mv .command.sh .${name}.command.sh
    mv .command.log .${name}.command.log
    """
}


process SYNTENY{
    label 'medium_resources'
    label 'syri'
    maxForks 10
    publishDir "${params.publishDir}/synteny", mode: 'copy'   

    input:
    tuple path(bam), path(bai)
    tuple path(reference), path(query)

    output:
    tuple path(".*command.sh"), path(".*command.log")

    script:
    """syri -c $bam -r $reference -q $query -k -F B --prefix ${bam.simpleName}
    plotsr ${bam.simpleName}syri.out $reference $query -H 8 -W 5
    """
//add markers and tracks --tracks tracks.txt --markers markers.bed
}

workflow {
    genome1 = Channel.fromPath(params.reference)
    genome2 = Channel.fromPath(params.query)
    chromosome_pairs = genome1.combine(genome2)
    alignments = MINIMAP2(chromosome_pairs)
    SYNTENY(alignments.mappings, chromosome_pairs)
    }
