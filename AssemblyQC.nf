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
    label 'kmertools'
    label 'medium_resources'
    publishDir "${params.publishDir}/assembly/QC/", mode: 'copy'
    errorStrategy 'ignore'

    input:
    tuple path(reads), path(genome)


    output:
    path "02_kat*"

    script:
    """kat comp $reads $genome -o kat -t $task.cpus
    mkdir 02_kat_${genome.baseName}
    cp .command.* 02_kat_${genome.baseName}
    mv kat* 02_kat_${genome.baseName}
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
    generate_plot.py -wd 01_busco_out_${genome.baseName}
    cp .command.* 01_busco_out_${genome.baseName} 
    """
}

process RUN_INSPECTOR {
    label 'inspector'
    label 'resource_intensive'
    publishDir "${params.publishDir}/assembly/QC/", mode: 'copy', pattern:"03_inspector_out*/"
   
    input:
    tuple path(reads), path(genome)


    output:
    path "03_inspector_out_*/"
    path '03_inspector_out_*/read_to_contig.bam.gz', emit: mappings

    
    script:
    """inspector.py -c $genome -r $reads -o 03_inspector_out_${genome.baseName} --datatype $params.dtype --thread $task.cpus --min_contig_length_assemblyerror 10000
    pigz 03_inspector_out_${genome.baseName}/read_to_contig.bam
    rm -r 03_inspector_out_${genome.baseName}/*_workspace
    cp .command.sh .command.log 03_inspector_out_${genome.baseName}/
    """
}

process DBLAST {
    label 'resource_intensive'
    label 'blobtools'

    input:
    path genome
    val db

    output
    path '*_uniref90.out', emit: blastout

    script:
    """
    cp $params.db uniref90.fasta.gz
    gunzip uniref90.fasta.gz
    diamond makedb --in $db --db uniref90
    diamond blastx --query $genome --db uniref90 \
    --outfmt 6 --sensitive --max-target-seqs 1 \
    --evalue 1e-25 --thread $task.cpus --out ${genome.baseName}_uniref90.out
    rm uniref90.fasta
    blobtools taxify -f ${genome.baseName}_uniref90.out -m $params.taxid -s 0 -t 1
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

process RETRIEVE_MITO {
    label 'low_resources'
    publishDir "${params.publishDir}/assembly/QC/05_mitochondria", mode: 'copy'

    input:
    path genome

    output:
    path '*_mito.out'
    path '.command.*'
    path 'putative_mito.fasta', emit: mitochondria

    script:
    """
    makeblastdb -dbtype nucl -in $genome
    blastn -query $projectDir'/lib/'$params.mito -db ${genome.name} -outfmt 6 -evalue 1e-25 -num_threads $task.cpus -max_target_seqs 5 -out ${genome.baseName}_mito.out
    sort -k3 -u *_mito.out | awk '{print \$2}' | head -1 > putative_contig
    seqkit grep -f putative_contig $genome > putative_mito.fasta"""
}

process CIRCULARIZE{
    publishDir "${params.publishDir}/assembly/QC/05_mitochondria", mode: 'copy'
    label 'nanoplot'
    label 'low_resources'
    input:
    path mitochondria

    output:
    path "circMitochondria.fasta"

    script:
    """
    simple_circularise.py $mitochondria circMitochondria.fasta -min 10000
    """
}

workflow CONTAMINATION {
    take:
        READS
        ASSEMBLIES
        MAPPINGS // puede estar vacio
    main:
        include { MINIMAP } from './modules/minimap.nf'
        
        if (mapped == false){
            MAPPINGS = MINIMAP(READS, ASSEMBLIES) // tuple output ASSEMBLY, BAM
            }
        else {
            MAPPINGS.Channel.fromPath(params.mapped)
            MAPPINGS.ifEmpty{MAPPINGS = MINIMAP(READS, ASSEMBLIES)}} 

        BLOBINPUT = DBLAST(MAPPINGS).blastout // tuple input ASSEMBLY, BAM // output tuple ASSEMBLY, BAM, BLASTOUT     
        BLOB = RUN_BLOBTOOLS(BLOBINPUT) // TUPLE INPUT; END
}

process MAKECFG {
    publishDir "${params.publishDir}/assembly/QC/05_mitochondria/mitonextpolish/", mode: 'copy', pattern: ".command*"
    
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
    publishDir "${params.publishDir}/assembly/QC/05_mitochondria", mode: 'copy', pattern: "*nextpolish"    

    input:
    path config_file
    path reads
    
    output:
    path 'mitonextpolish'
    path('mitonextpolish/mitochondria_polished.fasta'), emit: assembly

    script:
    """
    ls $reads > lgs.fofn
    nextPolish $config_file
    seqkit seq 01_rundir/genome.nextpolish.fasta -u > 01_rundir/mitochondria_polished.fasta
    mkdir mitonextpolish
    mv 01_rundir/mitochondria_polished.fasta mitonextpolish
    mv 01_rundir/*fasta.stat mitonextpolish
    cp .command.sh .command.log mitonextpolish
    """
}

workflow ASSEMBLE_MITOCHONDRIA{
    take:
        READS
        ASSEMBLIES
    main:
        MITO_BLAST = RETRIEVE_MITO(ASSEMBLIES).mitochondria
        config = MAKECFG(MITO_BLAST).config_file
        POLISHED = POLISH(config, READS).assembly
}

workflow ALL {
    take:
        READS
        ASSEMBLIES
        MAPPINGS
    main:
        INPUT = READS.combine(ASSEMBLIES)
        CONTINUITY = RUN_GFSTATS(INPUT)
        GCOMPLETENESS = RUN_BUSCO(INPUT)
        RCOMPLETENESS = RUN_KAT(INPUT)
        CORRECTNESS = RUN_INSPECTOR(INPUT)
        MAPPINGS = CORRECTNESS.mappings
            //CLEANESS = CONTAMINATION(READS, ASSEMBLIES, MAPPINGS)
}

workflow {
    READS = Channel.fromPath(params.reads, checkIfExists:true) 
    ASSEMBLIES = Channel.fromPath(params.genome, checkIfExists: true)    
    INPUT = READS.combine(ASSEMBLIES)
 
    MAPPINGS = channel.empty()    
    if (params.all == true){
        ALL(READS, ASSEMBLIES, MAPPINGS)}
    if (params.mitassembly == true) {
        ASSEMBLE_MITOCHONDRIA(READS, ASSEMBLIES)}
    if (params.continuity == true){
        CONTINUITY = RUN_GFSTATS(INPUT)}
    if (params.genecompleteness == true){
        GCOMPLETENESS = RUN_BUSCO(INPUT)}
    if (params.readcompleteness == true){
        RCOMPLETENESS = RUN_KAT(INPUT)}
    if (params.correctness == true){
        CORRECTNESS = RUN_INSPECTOR(INPUT)
        MAPPINGS = CORRECTNESS.mappings}
    if (params.cleaness ==true){
        CLEANESS = CONTAMINATION(READS, ASSEMBLIES, MAPPINGS)
    }
}

params.continuity = false
params.genecompleteness = false
params.readcompleteness = false
params.correctness = false
params.cleaness = false
params.all = false
params.mito ="Crenicichla_mitoGenes.fasta"
params.dtype = 'ont'
params.mitassembly = false
params.template2 = "/home/fhenning/assembly_project/lib/nextpolish.cfg"
params.taxid = "/home/fhenning/assembly_project/lib/uniref.taxid.tsv"
params.db = "/home/fhenning/assembly_project/lib/uniref90.fasta.gz"