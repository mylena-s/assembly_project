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
    publishDir "${params.publishDir}/annotation/", mode: 'copy'
    
    input:
    path genome

    output:
    path '01_earlGrey_*'

    script:
    """
    earlGrey -g $genome -s "" -o ./01_earlGrey_${genome.baseName} -t $task.cpus
    mv .command.log .command.sh 01_earlGrey_${genome.baseName}
    """
}    
params.proteins = "lib/Actinopterygii.fasta.gz" 

process BRAKER2 { 
    label 'braker'
    label 'resource_intensive'
    publishDir "${params.publishDir}/annotation/", mode: 'copy'
    
    input:
    path genome
    path proteins

    output:
    script:
    """
    gunzip -c $proteins > proteins.fa
    braker.pl --species=Crenicichla --genome=$genome --prot_seq=proteins.fa --workingdir=braker2_${genome.baseName} --threads=$task.cpus --busco_lineage=actinopterygii_odb10
    """
}


workflow {
    ASSEMBLIES = Channel.fromPath(params.genome, checkIfExists: true)
    REFPROT = Channel.fromPath(params.proteins,  checkIfExists: true)
    INPUT = ASSEMBLIES.combine(REFPROT)
    ANNOT = REPEAT_ANNOT(INPUT)
}