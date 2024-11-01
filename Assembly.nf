#!/usr/bin/env nextflow

log.info """\
    P I P E L I N E    F O R   A S S E M B L Y 
    ==========================================
    Uses NextDenovo as ONT assembler or Hifiasm
                 for Pacbio hifi. 
        
    Current project dir. is: $projectDir"""
    .stripIndent()

process PREPARE_CONFIG {
    label 'low_resources'
    
    input:
    val(template)

    output:
    path('config'), emit: config_file

    script:
    """#!/usr/bin/env python3
with open('$template', 'r') as file:
    content = file.read()
    content = content.replace(' = clr', ' = ont')
    content = content.replace('parallel_jobs = 2', 'parallel_jobs = 4')
    content = content.replace('genome_size = 308161', 'genome_size = ${params.gs}000000')
    content = content.replace('pa_correction = 2', 'pa_correction = 4')
    content = content.replace('correction_options = -p 15', 'correction_options = -p 12')
    content = content.replace('minimap2_options_raw =  -t 8', 'minimap2_options_raw =  -t 12')
with open('config', 'w') as file:
    file.write(content)
    """
}

process ASSEMBLY_ONT {
    label 'nextpolish'
    label 'resource_intensive'
    errorStrategy 'ignore'

    publishDir "${params.publishDir}/assembly/", mode: 'copy', pattern: "*nextdenovo"
    
    input:
    path reads
    path config

    output:
    path 'nextdenovo'
    path 'nextdenovo/nextdenovo.upper.fasta', emit: assembly

    script:
    """
    cp $config run.cfg
    ls $reads > input.fofn
    nextDenovo run.cfg
    seqkit seq 01_rundir/03.ctg_graph/*fasta -u > 01_rundir/03.ctg_graph/nextdenovo.upper.fasta
    mkdir nextdenovo
    mv 01_rundir/03.ctg_graph/*fasta nextdenovo
    mv 01_rundir/03.ctg_graph/*fasta.stat nextdenovo
    mv .command.sh .command.log nextdenovo
    cp run.cfg nextdenovo
    """
}

process MAKECFG {
    publishDir "${params.publishDir}/assembly/nextpolish/", mode: 'copy', pattern: ".command*"
    
    input:
    val(genome)
    
    output:
    path 'run.cfg', emit: config_file
    path '.command*'
    
    script:
    """#!/usr/bin/env python3
with open('$params.template2', 'r') as file:
  content = file.read()
  content = content.replace('INPUT','$genome')
with open('run.cfg', 'w') as file:
  file.write(content)
    """
}

process POLISH {
    label 'nextpolish'
    label 'medium_resources'
    publishDir "${params.publishDir}/assembly/", mode: 'copy', pattern: "*nextpolish"    

    input:
    path config_file
    path reads
    
    output:
    path 'nextpolish'
    path('nextpolish/nextpolish.upper.fasta'), emit: assembly

    script:
    """
    ls $reads > lgs.fofn
    nextPolish $config_file
    seqkit seq 01_rundir/genome.nextpolish.fasta -u > 01_rundir/nextpolish.upper.fasta
    mkdir nextpolish
    mv 01_rundir/*fasta nextpolish
    mv 01_rundir/*fasta.stat nextpolish
    cp .command.sh .command.log nextpolish
    """
}

process  ASSEMBLY_HIFI {
    label 'resource_intensive'
    publishDir "${params.publishDir}/assembly/hifiasm/", mode: 'copy', pattern: "hifiasm/"

    input:
    path reads
    path ul 

    output:
    path "hifiasm"
    path "hifiasm/*asm.bp.p_ctg.gfa.fasta", emit: assembly
    
    script:
    def filter = ul != 'NO_FILE' ? "--ul $ul" : ''
    """
    hifiasm -o ${reads.baseName}.asm -t $task.cpus $filter $reads
    mkdir hifiasm 
    mv ${reads.baseName}.asm* hifiasm
    mv .command* hifiasm
    cd hifiasm
    for file in *gfa; do gfastats $file -o fa --discover-paths > $file'.fasta'
    """
}

process ASSEMBLY_HIBRID{
    label 'verkko'
    label 'resource_intensive'
    publishDir "${params.publishDir}/assembly/", mode: 'copy'

    input:
    tuple path(ONT_reads), path(HIFI_reads)
    
    output:
    path "verkko", emit: folder
    path "verkko/"
    
    script:
    """
    verkko -d ./verkko --hifi $HIFI_reads --nano $ONT_reads --local-cpus $task.cpus
    cp .command* verkko
    """
}

process BREAK_MISSASSEM {
    label 'goldrush'
    label 'medium_resources'

    publishDir "${params.publishDir}/assembly/tigmint/", mode: 'copy'

    input:
    path genome, name: "genome.fa"
    path reads, name: "reads.fq.gz"

    output:
    path "breaktigs.fa", emit: output
    path "*.bed"
    path "*.breaktigs.fa"
    path "nametable"
    script:
    
    """tigmint-make tigmint-long draft=genome reads=reads span=auto G=${params.gs}000000 dist=auto
    grep ">" *.breaktigs.fa | sed 's/>//g' - > nameTable.temp 
    seqkit replace -p '.+' -r 'ctg_{nr}' *.breaktigs.fa > breaktigs.fa 
    grep ">" breaktigs.fa | sed 's/>//g' - > nameTable2.temp
    paste -d '\t' nameTable.temp nameTable2.temp > nametable 
    """ 
}

process SCAFFOLDING {
    label 'goldrush'
    label 'medium_resources'
    publishDir "${params.publishDir}/assembly/ntlink_scaffolding/", mode: 'copy'  
    errorStrategy 'ignore'

    input:
    path genome 
    path reads,  name: "reads.fq.gz"

    output:
    path("ntlink.fasta"), emit: scaffolded_genome
    path "ntlink.agp", emit: assemblygraph
    path "*paf", emit: alignment
    path ".command*"

    script:
    """
    seqkit seq -m 1000 $genome > genome_subset.fasta
    ntLink_rounds run_rounds_gaps target=genome_subset.fasta reads=reads.fq.gz k=32 w=250 t=$task.cpus rounds=4 clean=True overlap=True a=$params.nreads paf=True
    mv genome_subset.fasta.k32.w250.z1000.ntLink.ntLink.gap_fill.fa.k32.w250.z1000.ntLink.scaffolds.gap_fill.fa ntlink.fasta
    mv genome_subset.fasta.k32.w250.z1000.ntLink.ntLink.gap_fill.fa.k32.w250.z1000.ntLink.scaffolds.gap_fill.fa.agp ntlink.agp
    """
}


process PURGE_HAPLOTIGS {
    label 'haplotigs'
    label 'medium_resources'
    publishDir "${params.publishDir}/assembly/", mode: 'copy'

    input:
    path genome
    path reads
    val type

    output:
    path "haplotigs_out_${val}"
    path "haplotigs_out_${val}/curated.fasta", emit: curated

    script:
    """
    minimap2 -ax map-$type $genome $reads --secondary=no | samtools sort -o mapping_LR.map-${type}.bam -T tmp.ali
    purge_haplotigs hist -b mapping_LR.map-${type}.bam  -g $genome  -t $task.cpus
    purge_haplotigs cov -i mapping_LR.map-${type}.bam.200.gencov -l 5 -m 190 -h 190 -o coverage_stats.csv -j 190 -s 80
    purge_haplotigs purge -g $genome -c coverage_stats.csv -t $task.cpus -d -b mapping_LR.map-${type}.bam
    mkdir haplotigs_out_${genome.baseName}
    mv curated* coverage_stats.csv dotplots* mapping_LR.map-${type}.bam* haplotigs_out_${genome.baseName}
    cp .command.sh .command.log haplotigs_out_${genome.baseName}"""
}

params.template1 = "/home/fhenning/assembly_project/lib/nextdenovo.cfg"
params.template2 = "/home/fhenning/assembly_project/lib/nextpolish.cfg"
params.assembled = false 
params.type = "ont" 
params.gs = 900
params.db = '/home/fhenning/assembly_project/lib/actinopterygii_odb10.tar.gz'
params.nreads = 2
params.ul = "NO_FILE"


include { SAMPLE_ONT } from './modules/extra.nf'

workflow {
    READS = Channel.fromPath(params.reads, checkIfExists:true) 
    UL = Channel.fromPath(params.ul, checkIfExists:false)
    if (params.assembled == false) {
        if (params.type == "ont") {
            CONFIG = PREPARE_CONFIG(params.template1)
            PRIMARY = ASSEMBLY_ONT(READS,CONFIG).assembly
            BREAKTIGS = BREAK_MISSASSEM(PRIMARY, READS).output
            HAPLOTIGS = PURGE_HAPLOTIGS(BREAKTIGS, READS, "ont").curated
            SCAFFOLD = SCAFFOLDING(HAPLOTIGS, READS).scaffolded_genome
            config = MAKECFG(SCAFFOLD).config_file
            POLISHED = POLISH(config, READS).assembly
            ASSEMBLIES = PRIMARY.concat(BREAKTIGS, HAPLOTIGS, SCAFFOLD, POLISHED)

        } else {
            if (params.type == "hifi"){
                PRIMARY = ASSEMBLY_HIFI(READS, UL).assembly
                HAPLOTIGS = PURGE_HAPLOTIGS(PRIMARY, READS, params.type)

            } else {
                ONT=SAMPLE_ONT(UL, 50000).reads
                INPUT = ONT.combine(READS)
                PRIMARY = ASSEMBLY_HIBRID(INPUT)
            }
    }
    }
    else {
        PRIMARY = Channel.fromPath(params.assembled)
        HAPLOTIGS = PURGE_HAPLOTIGS(PRIMARY, READS, params.type)

    }
}