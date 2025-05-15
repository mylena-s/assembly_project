process MAP{
    label 'medium_resources'
    label 'methylation'
    maxForks 2
    publishDir "${params.publishDir}/methylation", mode: 'copy'   

    input:
    tuple path(genome), path(bam)
    
    output:
    path("methylation.bam"), emit: bam

    script:
"""gunzip -c $bam > reads.bam 
pbmm2 align $genome reads.bam methylation.bam
mv .command.sh .methylation.command.sh
mv .command.log .methylation.command.log
"""
}


process CALL{
    label 'medium_resources'
    label 'syri'
    maxForks 10
    publishDir "${params.publishDir}/methylation", mode: 'copy'   

    input:
    path(bam)
    output:
    path("${bam.simpleName}*")
    path("sorted.bam")

    script:
    """
samtools sort $bam > sorted.bam
samtools index sorted.bam
aligned_bam_to_cpg_scores --bam sorted.bam --output-prefix ${bam.simpleName} --model $projectDir'/bin/pileup_calling_model.v1.tflite' --threads $task.cpus
"""
}



params.aligned = false
workflow {
    genome = Channel.fromPath(params.reference)
    if(params.aligned == false){
        reads_bam = Channel.fromPath(params.bam)
        input = genome.combine(reads_bam)
        alignments = MAP(input).bam
    } else {
        alignments = Channel.fromPath(params.bam)}    
    calls = CALL(alignments)
    }
