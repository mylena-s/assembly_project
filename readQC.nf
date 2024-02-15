#!/usr/bin/env nextflow

log.info """\
    P I P E L I N E    F O R   A S S E M B L Y 
    ==========================================
    QC step for reads

    Uses NanoPlot to assess read quality
    porechop and cutadapt for trimming (ont and hifi)
    and Genomescope for pre-assembly assessments
    (genome size, heterozigosity etc)"""
    .stripIndent()


process QC {
    publishDir "${params.publishDir}/pre_assembly/nanoplot/", mode: 'copy'
    label 'nanoplot'
    label 'low_resources'
    
    input:
    path reads

    output:
    path 'NanoPlot*', emit: folder
  
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
    path "*_trimmed.fq.gz", emit: reads
    path ".command.sh"
    path ".command.log", emit: report

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

process CUTADAPT{
    label 'cutadapt'
    label 'low_resources'

    input:
    path reads
    output:
    path '*_trimmed.fq.gz', emit: reads
    path '*.cutadapt.json', emit: report
    script:
    """
    cutadapt -b ATCTCTCTCAACAACAACAACGGAGGAGGAGGAAAAGAGAGAGAT \
    -b ATCTCTCTCTTTTCCTCCTCCTCCGTTGTTGTTGTTGAGAGAGAT \
    -e 0.1 -O 35 --rc --discard-trimmed -j 0 \
    -o ${reads.baseName}_trimmed.fq.gz $reads \
    -j 10 --json=${reads.baseName}.cutadapt.json
    """
}

process MULTIQC{
    label 'multiqc'
    input:
    path reports
    output:
    path 'multiqc'
    script:
    """
    multiqc .
    """
}
workflow {
    SCRUBBED = channel.empty()    
    REPORTS =  channel.empty()
    if (params.trimmed == false) {
        READS = Channel.fromPath(params.reads, checkIfExists:true)
        if (params.type == 'ont') {
            RESULTS = ADAPTOR_CHECK(READS)
            TRIMMED = RESULTS.out.reads
            REPORTS = RESULTS.report.out.collect()}
            
        if (params.type == 'hifi') {
            RESULTS = CUTADAPT(READS)}
            TRIMMED = RESULTS.reads
            REPORTS = RESULTS.report.collect()
    } else {
        READS = channel.empty()    
        TRIMMED = Channel.fromPath(params.reads, checkIfExists:true)}
    if (params.scrub == true) {
        SCRUBBED = CHIMERA_CHECK(TRIMMED).reads}

    READS_CH = READS.concat(READS, TRIMMED, SCRUBBED)
    NANOPLOT = QC(READS_CH).folder.collect()

    if (params.statistics == true){
        KMER_HIST = RUN_JELLYFISH(READS)
        REPORT = GENOMESCOPE(KMER_HIST.hist).genomescope_report}
}

params.type = "ont"
params.kmer = 21
params.trimmed = false
params.scrub = false
params.statistics = false
