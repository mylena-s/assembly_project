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
    publishDir "${params.publishDir}/assembly/QC/", mode: 'copy', saveAs:"03_inspector_out_${genome.baseName}", pattern:"03_inspector_out/"
   
    input:
    tuple path(reads), path(genome)


    output:
    path '03_inspector_out'
    path '03_inspector_out/read_to_contig.bam', emit: mappings

    
    script:
    """inspector.py -c $genome -r $reads -o 03_inspector_out --datatype $params.type --thread $task.cpus --min_contig_length_assemblyerror 10000
    pigz 03_inspector_out/read_to_contig.bam
    rm -r 03_inspector_out/*_workspace
    cp .command.sh .command.log 03_inspector_out
    """
}

process BLAST {
    label 'resource_intensive'
    
    input:
    path genome

    val db

    output
    path ' blobtools_blast.out', emit: blastout

    script:
    """
    blastn -query $genome -db $db -outfmt "6 qseqid staxids bitscore std" -num_threads $task.cpus -evalue 1e-25 \
    -max_target_seqs 10 -max_hsps 1 -out blobtools_blast.out
"""
}

process MINIMAP{
    label 'medium_resources'
    label 'inspector'

    input:
    path genome
    path reads

    output:
    path 'read_to_contig.bam', emit: mappings

    script:
    """minimap2 -ax map-$params.dtype $genome $reads | samtools sort -@$task.cpus -O BAM -o read_to_contig.bam -
    """
}

process RUN_BLOBTOOLS {
    label 'medium_resources'
    label 'blobtools'

    input:
    path genome
    path blastout
    path sortedbam
    
    output:
    path 'blobtools_dir', emit: folder
    tuple path('blobtools_dir/blobDB.table.txt'), val ("${transcriptome.baseName}"), emit: blobout
   
    script:
    """
    blobtools create --fasta $genome 04_blobtools_${genome.baseName}
    blobtools add --hits $blastout --cov $sortedbam --taxdump $projectDir/lib/taxdump/ 04_blobtools_${genome.baseName}
    """    
}
params.continuity = true
params.genecompleteness = false
params.readcompleteness = false
params.correctness = false
params.cleaness = false
params.all = false
params.dtype = 'ont'

workflow CONTAMINATION {
    take:
        READS
        ASSEMBLIES
        MAPPINGS // puede estar vacio
    main:
        if (mapped == false){
            MAPPINGS = MINIMAP(READS, ASSEMBLIES) // tuple output ASSEMBLY, BAM
            }
        else {
            MAPPINGS \
            | ifEmpty { 'EMPTY' } \
            | MAPPINGS = MINIMAP(READS, ASSEMBLIES)}
        BLOBINPUT = BLAST(MAPPINGS).blastout // tuple input ASSEMBLY, BAM // output tuple ASSEMBLY, BAM, BLASTOUT     
        BLOB = RUN_BLOBTOOLS(BLOBINPUT) // TUPLE INPUT; END
}

workflow {
    ASSEMBLIES = Channel.fromPath(params.genome, checkIfExists: true)    
    INPUT = READS.combine(ASSEMBLIES)
    MAPPINGS = channel.empty()    

    if (params.continuity == true){
        CONTINUITY = RUN_GFSTATS(INPUT)}
    if (params.genecompleteness == true){
        GCOMPLETENESS = RUN_BUSCO(INPUT)}
    if (params.readcompleteness == true){
        RCOMPLETENESS = RUN_KAT(INPUT)}
    if (params.correctness == true){
        CORRECTNESS = RUN_INSPECTOR(INPUT)}
        MAPPINGS = CORRECTNESS.mappings
//    if (params.cleaness == true){
//     CLEANESS = CONTAMINATION()}

    if (params.all == true){
                CONTINUITY = RUN_GFSTATS(INPUT)
                GCOMPLETENESS = RUN_BUSCO(INPUT)
                RCOMPLETENESS = RUN_KAT(INPUT)
                CORRECTNESS = RUN_INSPECTOR(INPUT)
                MAPPINGS = CORRECTNESS.mappings  ///// tuple output ASSEMBLY, BAM 
                CLEANESS = CONTAMINATION() 
}
}