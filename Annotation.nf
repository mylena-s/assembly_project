#!/usr/bin/env nextflow

log.info """\
    P I P E L I N E    F O R   R E P E A T 
              A N N O T A T I O N 
    =======================================
    Uses NextDenovo as ONT assembler or Hifiasm
                 for Pacbio hifi. 
    
    Also performs QC with different tools
    
    Current project dir. is: $projectDir"""
    .stripIndent()

process REPEAT_ANNOT {
    label 'earlGrey'
    label 'resource_intensive'
    publishDir "${params.publishDir}/annotation/", mode: 'copy', pattern: "*01_earlGrey"
    
    input:
    path genome

    output:
    path '01_earlGrey*'

    script:
    """
    earlGrey -g $genome -s "" -o ./01_earlGrey -t $task.cpus
    mv .command.log .command.sh 01_earlGrey
    """
}    

params.gs = 900

workflow {
    ASSEMBLIES = Channel.fromPath(params.genome, checkIfExists: true)
    ANNOT = REPEAT_ANNOT(ASSEMBLIES)
}