#!/usr/bin/env nextflow

log.info """\
    P I P E L I N E    F O R   A S S E M B L Y 
    ==========================================
    QC step for reads

    Uses NanoPlot to assess read quality
    and Genomescope for pre-assembly assessments
    (genome size, heterozigosity etc)"""
    .stripIndent()


process QC {
    publishDir "${params.publishDir}/pre_assembly/nanoplot/", mode: 'copy'
    label 'nanoplot'
    label 'medium_resources'
    
    input:
    path reads

    output:
    path 'NanoPlot*', emit: output
  
    script:
    """
    NanoPlot -t $task.cpus --verbose --huge -o NanoPlot_${reads.baseName} --tsv_stats --raw --fastq $reads
    """
}

process GENOMESCOPE {
    publishDir "${params.publishDir}/pre_assembly/genomescope/", mode: 'copy'    
    label "low_resources"
    label "genomescope"

    input:
    path hist

    output:
    path 'report',      emit: genomescope_report
  
    script:
    """
    genomescope.R -i $hist -o report -k $params.kmer 
    """
}

process RUN_JELLYFISH{
    publishDir "${params.publishDir}/pre_assembly/genomescope/", mode: 'copy'    
    label 'low_resources'

    input:
    path reads
    
    output:
    path '*.histo', emit:hist
    
    script:
    """
    gunzip -c $reads > reads.fastq 
    jellyfish count -C -m $params.kmer -s 1000000000 -t $task.cpus reads.fastq -o reads.jf
    rm reads.fastq
    jellyfish histo -t $task.cpus reads.jf > ${reads.baseName}.histo
    rm reads.jf """
}

process ADAPTOR_CHECK {
    label 'porechop'
    label 'medium_resources'
    errorStrategy 'ignore'
    publishDir "${params.publishDir}/pre_assembly/trimming/", mode: 'copy'    


    input:
    path reads

    output:
    path "*_trimmed.fq.gz", emit: reads, optional: true
    path ".command.sh"
    path ".command.log"

    script:
    """
    porechop_abi -abi -i $reads -tmp ./temp -o ${reads.baseName}_trimmed.fq -t $task.cpus
    pigz *_trimmed.fq"""
}

process CHIMERA_CHECK {
    label 'yacrd'
    label 'medium_resources'
    errorStrategy 'ignore'
    publishDir "${params.publishDir}/pre_assembly/chimera/", mode: 'copy'    

    input:
    path reads
    output:
    path "*_scrubb.fastq.gz", emit: reads
    path ".command.sh"
    path ".command.log"

    script:
//thre is an error in minimap command
    """minimap2 -x ava-ont -g 500 $reads $reads > mapping.paf
    yacrd -i overlap.paf -o reads.yacrd
    yacrd -i mapping.paf -o reads.yacrd scrubb -i $reads -o ${reads.baseName}_scrubb.fastq.gz
    """
}

workflow {
    if (params.trimmed == false) {
        READS = Channel.fromPath(params.reads, checkIfExists:true)
        TRIMMED = ADAPTOR_CHECK(READS).reads
    } else {
       READS = channel.empty()    
       TRIMMED = Channel.fromPath(params.reads, checkIfExists:true)
    }

    SCRUBBED = CHIMERA_CHECK(TRIMMED).reads
    READS_CH = READS.concat(TRIMMED, SCRUBBED)
    QC(READS_CH)
//    KMER_HIST = RUN_JELLYFISH(READS)
//    REPORT = GENOMESCOPE(KMER_HIST.hist).genomescope_report
}

params.type = "ONT"
params.kmer = 15
params.trimmed = false
