
log.info """\
    P I P E L I N E     F O R    P O P G E N  
    ==========================================
   
    performs long read mapping with minimap
    sorts, marks duplicates, and calls variants    
    Current project dir. is: $projectDir"""
    .stripIndent()

include { MINIMAP; BAMQC; FREEBAYES_JOINT; MARKDUP } from './modules/minimap.nf'  
params.dtype = "hifi"
params.mapped = false
params.nodup = false

process JOIN_VCF {
    label 'low_resources'
    label 'minimap'
    publishDir "${params.publishDir}/population/SNP_variant_calling", mode: 'copy'   

    input:
    path(vcfs)

    output:
    path "all_jointcall.vcf.gz"
    path "concat.command.*"

    script:
    def files = vcfs.join(' ')
    """
    bcftools concat $files > all_jointcall.vcf
    bgzip all_jointcall.vcf
    mv .command.sh concat.command.sh
    mv .command.log concat.command.log
    """
}

process MITOHIFI{
    label 'medium_resources'
    label 'mitohifi'
    errorStrategy 'ignore'

    publishDir "${params.publishDir}/population/mitochondrial", mode: 'copy'   

    input:
    path read
    val reference

    output:
    path "*final_mitogenome.*"
    path ".*command.*"
    
    script:
    """
    mitohifi.py -r $read -f ${reference}.fasta -g ${reference}.gb -t $task.cpus -a animal 
    mv final_mitogenome.fasta ${reference}.final_mitogenome.fasta
    mv final_mitogenome.gb ${reference}.final_mitogenome.gb
    mv .command.log .${reference}.command.log
    mv .command.sh .${reference}.command.sh
    """
}

reference = "$projectDir/lib/NC_030272.1"
params.assembled_mito = false

workflow MITO {
    take:
    READS
    main:
    if (params.assembled_mito == false){
        MITOCONDRIAS = MITOHIFI(READS,reference)}
    else {
        MITOCONDRIAS = Channel.fromPath(params.assembled_mito)}
    // ALIGN
    // 
    }


workflow NUCLEAR {    
    take:
    READS

    main:
    ASSEMBLY = Channel.fromPath(params.genome, checkIfExists: true)    
    INPUT = ASSEMBLY.combine(READS)
    COLLECTED_MAPPINGS = channel.empty()

    if (params.mapped == false) {
        MAPPINGS = MINIMAP(INPUT).mappings
    } else {
        MAPPINGS = Channel.fromPath(params.mapped, checkIfExists: true) 
    }
    if (params.nodup == false) {
        MAPPINGS_NODUP = MARKDUP(MAPPINGS).mappings
        COLLECTED_MAPPINGS = MAPPINGS_NODUP.collect()
    } else {
        MAPPINGS_NODUP = Channel.fromPath(params.nodup, checkIfExists:true) }
        COLLECTED_MAPPINGS = MAPPINGS_NODUP.collect() 
//    BAMQC(MAPPINGS_NODUP)
    contigs_file = file(params.contigs)
    allcontigs = Channel.fromList(contigs_file.readLines())
    GENOME_CONTIG = ASSEMBLY.combine(allcontigs)
    if (params.vcfs == false) {
        VCFS = FREEBAYES_JOINT(GENOME_CONTIG, COLLECTED_MAPPINGS).vcf.collect()
    } else {
        VCFS = Channel.fromPath(params.vcfs, checkIfExists: true).collect()}   
    ALL_VCF = JOIN_VCF(VCFS)
}


workflow {
    READS = Channel.fromPath(params.reads, checkIfExists:true)
    // NUCLEAR(READS)
    MITO(READS)
    }
    
    
