#!/usr/bin/env nextflow

log.info """\
    P I P E L I N E    F O R   A S S E M B L Y 
    ==========================================
    Uses NextDenovo as ONT assembler or Hifiasm
                 for Pacbio hifi. 
    
    Also performs QC with different tools
    
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
    tuple val(genome)
    
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
    script:
    ""
}

process SCAFFOLDING {
    label 'goldrush'
    label 'medium_resources'
    publishDir "${params.publishDir}/assembly/ntlink_scaffolding/", mode: 'copy'  
    errorStrategy 'ignore'

    input:
    path genome
    path reads

    output:
    path("ntlink.fasta"), emit: scaffolded_genome
    path "ntlink.agp", emit: assemblygraph
    path ".command*"

    script:
    """
    ntLink_rounds run_rounds_gaps target=$genome reads=$reads k=32 w=250 t=$task.cpus rounds=3 clean=True overlap=True a=$params.nreads
    mv *.fasta.k32.w250.z1000.ntLink.ntLink.gap_fill.fa.k32.w250.z1000.ntLink.scaffolds.gap_fill.fa ntlink.fasta
    mv *.fasta.k32.w250.z1000.ntLink.ntLink.gap_fill.fa.k32.w250.z1000.ntLink.scaffolds.gap_fill.fa.agp ntlink.agp
    """
}


process PURGE_HAPLOTIGS {
    label 'haplotigs'
    label 'medium_resources'
    publishDir "${params.publishDir}/assembly/", mode: 'copy'    

    input:
    path genome
    path reads

    output:
    path "haplotigs_out"

    script:
    """
    minimap2 -ax map-ont $genome $reads --secondary=no | samtools sort -o mapping_LR.map-ont.bam -T tmp.ali
    purge_haplotigs hist -b mapping_LR.map-ont.bam  -g $genome  -t $task.cpus
    purge_haplotigs cov -i mapping_LR.map-ont.bam.gencov -l 5 -m 190 -h 190 -o coverage_stats.csv -j 190 -s 80
    purge_haplotigs purge -g $genome -c coverage_stats.csv -t $task.cpus
    mkdir haplotigs_out
    mv curated* coverage_stats.csv dotplots* mapping_LR.map-ont.bam.* haplotigs_out
    cp .command.sh .command.log haplotigs_out"""
}

///// STATS ///////

process RUN_GFSTATS {
    label 'low_resources'
    publishDir "${params.publishDir}/assembly/QC/00_gfstats_${genome.baseName}/", mode: 'copy'
    publishDir "${params.publishDir}/assembly/QC/00_gfstats_${genome.baseName}/", mode: 'copy'
    publishDir "${params.publishDir}/assembly/QC/00_gfstats_${genome.baseName}/", mode: 'copy'
    
    input:
    tuple path(reads), path(genome)

    output:
    path '*.stats'
    path '.command.*'

    script:
    gs = params.gs * 100000
    """
    gfastats $genome $params.gs --nstar-report -t > ${genome.baseName}.stats
    """
}

process RUN_KAT {
    label 'kat'
    label 'medium_resources'
    publishDir "${params.publishDir}/assembly/QC/", mode: 'copy', pattern:"02_kat", saveAs: "02_kat_${genome.baseName}"    
    errorStrategy 'ignore'

    input:
    tuple path(reads), path(genome)


    output:
    path "02_kat*"

    script:
    """kat comp $reads $genome -o kat -t $task.cpus
    mkdir 02_kat
    mv .command.* 02_kat
    mv kat* 02_kat
    cp .command.sh .command.log 02_kat
    """
}

process RUN_BUSCO {
    label 'busco'
    label 'medium_resources'
    publishDir "${params.publishDir}/assembly/QC/", mode: 'copy'
    
    input:
    tuple path(reads), path(genome)


    output:
    path '01_busco_out*'

    script:

    """
    busco -i $genome -l actinopterygii_odb10 -o 01_busco_out_${genome.baseName} -m genome -c $task.cpus
    cp .command.* 01_busco_out_${genome.baseName} 
    """
}

process RUN_INSPECTOR {
    label 'inspector'
    label 'resource_intensive'
    publishDir "${params.publishDir}/assembly/QC/", mode: 'copy', saveAs:"03_inspector_out_${genome.baseName}"    
   
    input:
    tuple path(reads), path(genome)


    output:
    path '03_inspector_out'

    
    script:
    """inspector.py -c $genome -r $reads -o 03_inspector_out --datatype $params.type --thread $task.cpus --min_contig_length_assemblyerror 10000
    gzip 03_inspector_out/read_to_contig.bam
    rm -r 03_inspector_out/*_workspace
    cp .command.sh .command.log 03_inspector_out
    """
}

process PREPARE_BLOBTOOLS {
    script:
    """blastn -task megablast -query $genome -db $database -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -culling_limit 5 -num_threads 50 -evalue 1e-25 -max_target_seqs 5 -out blobtools_blast.out
    minimap2 -ax map-ont ref.fa ont.fq.gz > aln.sam         # Oxford Nanopore genomic reads

"""
}


process RUN_BLOBTOOLS {
    label 'medium_resources'
    label 'blobtools'

    input:
    path genome
    path blastout
    tuple path(sortedbam), path(bamindex)
    
    output:
    path 'blobtools_dir', emit: folder
    tuple path('blobtools_dir/blobDB.table.txt'), val ("${transcriptome.baseName}"), emit: blobout
   
    script:
    """
    mkdir blobtools_dir
    blobtools create -i $transcriptome -t $blastout -b $sortedbam
    blobtools plot -i blobDB.json -r superkingdom
    blobtools plot -i blobDB.json -r phylum
    blobtools plot -i blobDB.json -r order
    blobtools view -i blobDB.json -r superkingdom
    mv *.png *txt blobtools_dir"""
}

workflow QC {
    take:
        INPUT

    main:
        CONTINUITY = RUN_GFSTATS(INPUT)
        COMPLETENESS = RUN_BUSCO(INPUT)
        COMPLETENESS = RUN_KAT(INPUT)
//        CORRECTNESS = RUN_INSPECTOR(INPUT)
//     CLEANESS = RUN_BLOBTOOLS()
}

params.template1 = "/home/fhenning/Mylena/lib/nextdenovo.cfg"
params.template2 = "/home/fhenning/Mylena/lib/nextpolish.cfg"
params.assemble = false 
params.type = "nanopore" 
params.gs = 900
params.db = '/home/fhenning/Mylena/lib/actinopterygii_odb10.tar.gz'
params.nreads = 5



workflow {
    READS = Channel.fromPath(params.reads, checkIfExists:true) 
    if (params.assemble != false) {
        if (params.type == "nanopore") {
            CONFIG = PREPARE_CONFIG(params.template1)
            PRIMARY = ASSEMBLY_ONT(READS,CONFIG).assembly
            SCAFFOLD = SCAFFOLDING(PRIMARY, READS).scaffolded_genome
            config = MAKECFG(SCAFFOLD).config_file
            POLISHED = POLISH(config, READS).assembly
            HAPLOTIGS = PURGE_HAPLOTIGS(POLISHED, READS)
            ASSEMBLIES = PRIMARY.concat(SCAFFOLD, POLISHED, HAPLOTIGS)

        } else {
            PRIMARY = ASSEMBLY_HIFI()
        }
    }
    // or read assembly/ies to do QC
    if (params.assemble == false) {
        POLISHED = Channel.fromPath(params.genome, checkIfExists: true)
        HAPLOTIGS = PURGE_HAPLOTIGS(POLISHED, READS)
        ASSEMBLIES = POLISHED.concat(HAPLOTIGS)
    }
    INPUT = READS.combine(ASSEMBLIES)
    //QC(INPUT) // MODIFY: QC SHOULD WORK WITH MORE THAN ONE ASSEMBLY

}