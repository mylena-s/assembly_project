#!/usr/bin/env nextflow

log.info """\
    P I P E L I N E    F O R   A S S E M B L Y 
    ==========================================
   
    performs assembly QC with different tools
    
    Current project dir. is: $projectDir"""
    .stripIndent()

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
    """
    gfastats $genome ${params.gs}00000 --nstar-report -t > ${genome.baseName}.stats
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
params.continuity = true
params.genecompleteness = false
params.readcompleteness = false
params.correctness = false
params.cleaness = false

workflow {
    ASSEMBLIES = Channel.fromPath(params.genome, checkIfExists: true)    
    INPUT = READS.combine(ASSEMBLIES)
    if (params.continuity == true){
        CONTINUITY = RUN_GFSTATS(INPUT)}
    if (params.genecompleteness == true){
        COMPLETENESS = RUN_BUSCO(INPUT)}
    if (params.readcompleteness == true){
        COMPLETENESS = RUN_KAT(INPUT)}
    if (params.correctness == true){
        CORRECTNESS = RUN_INSPECTOR(INPUT)}
//    if (params.cleaness == true){
//     CLEANESS = RUN_BLOBTOOLS()}
}