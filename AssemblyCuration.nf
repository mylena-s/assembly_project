#!/usr/bin/env nextflow

log.info """\
    P I P E L I N E    F O R   A S S E M B L Y 
    ==========================================
    Uses NextDenovo as ONT assembler or Hifiasm
                 for Pacbio hifi. 
        
    Current project dir. is: $projectDir"""
    .stripIndent()

process SEQTK {
    label 'low_resources'
    label 'kmertools'
    publishDir "${params.publishDir}/assembly", mode: 'copy'   

    input:
    path(genome)
    each path(list)

    output:
    path "${list.baseName}.fasta", emit: genome_subset

    script:
    """
    seqtk subseq $genome $list > ${list.baseName}.fasta
    """
}

include { SAMPLE_ONT } from './modules/extra.nf'
include { MINIMAP2 as MINIMAP} from './modules/minimap.nf'
include { MINIMAP2 as MINIMAP2 } from './modules/minimap.nf'

include { SCAFFOLDING } from './modules/extra.nf'
params.nreads = 1

workflow {
    GENOME = Channel.fromPath(params.genome, checkIfExists:true) 
    UL = Channel.fromPath(params.ul, checkIfExists:false)
    ONT = SAMPLE_ONT(UL, 100000)
    LISTS = Channel.fromPath(params.chr_list, checkIfExists:true)
    SAMPLED_GENOME = SEQTK(GENOME, LISTS).genome_subset
    SCAFFOLD = SCAFFOLDING(SAMPLED_GENOME, ONT)
    SCAFFOLDED_GENOME = SCAFFOLD.scaffolded_genome
    SCAFFOLDvsREADS = SCAFFOLDED_GENOME.combine(ONT)
    MAPPING = MINIMAP2(SCAFFOLDvsREADS, "ont").mappings
    SCAFFOLDvsCONTIGS = SCAFFOLDED_GENOME.combine(SAMPLED_GENOME)
    MINIMAP(SCAFFOLDvsCONTIGS, "hifi")
}