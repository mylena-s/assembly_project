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
    publishDir "${params.publishDir}/pre_assembly/", mode: 'copy'    

    input:
    path reads

    output:
    path 'NanoPlot',      emit: output
  
    script:
    """
    NanoPlot -t$task.cpus --verbose --huge -o NanoPlot --tsv_stats --raw --fastq $reads
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

workflow {
    READS = Channel.fromPath(params.reads, checkIfExists:true) 
//    QC()
    KMER_HIST = RUN_JELLYFISH(READS)
    REPORT = GENOMESCOPE(KMER_HIST.hist).genomescope_report
}

params.publishDir = ""
params.reads = ""
params.mem = 100
params.type = "ONT"
params.kmer = 15