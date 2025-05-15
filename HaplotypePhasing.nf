process PHASE_REFERENCE{
    label 'whatshap'
    label "resource_intensive"
    publishDir "${params.publishDir}/phasing/", mode: 'copy'

    input:
    path (reference)
    path(vcf)
    path(bam)

    output:
    path("phased.vcf.gz"), emit: phased_vcf
    path("phased.gtf"), emit: blocks
    path("*phasing_statistics.tsv")
    tuple path(".whatshap.command.sh"), path(".whatshap.command.log")

    script:
    """
    samtools faidx $reference
    samtools index -M $bam -@$task.cpus
    whatshap phase -o phased.vcf.gz --ignore-read-groups --reference=$reference $vcf $bam
    whatshap stats phased.vcf.gz --tsv phasing_statistics.tsv --gtf=phased.gtf
    tabix phased.vcf.gz
    bcftools consensus -H 1 -f $reference phased.vcf.gz > haplotype1.fasta
    bcftools consensus -H 2 -f $reference phased.vcf.gz > haplotype2.fasta                                                               
    mv .command.sh .whatshap.command.sh
    mv .command.log .whatshap.command.log
    """
}

process PANGENOME{
    label 'cactus'
    publishDir "${params.publishDir}/phasing/", mode: 'copy'
    label "resource_intensive"

    input:
    path genomes
    path reference

    output:
    path "pangenome"

    script:
    """cactus-pangenome $genomes --outDir pangenome --outName ${reference.baseName}_pangenome --reference ${reference.baseName} --vcf --giraffe --gfa --gbz"""
}

process REFERENCE_VCF{
    label 'medium_resources'
    publishDir "${params.publishDir}/phasing", mode: 'copy'   
    label 'vcflib'
    maxForks 10

    input:
    path vcf
    tuple val(mindep), val(maxdep)

    output:
    path "*.vcf.gz", emit: vcf
    tuple path(".command.sh"), path(".command.log")

    script:
    def name = "${vcf}".replaceAll(".vcf.gz","")
    """
    vcffilter -f 'QUAL / AO > 10' -g 'DP > $mindep & DP < $maxdep' $vcf > filtered.vcf
    bgzip filtered.vcf
    """
}

include { MINIMAP2; FREEBAYES_SINGLE } from './minimap.nf'
include { MINIMAP2 as MAP_ONT } from './minimap.nf'

params.aligned = false
params.ont = false
params.variant_call = true
params.filter = true
workflow {
    genome = Channel.fromPath(params.reference)
    if (params.variant_call != true){
            if(params.aligned == false){
                hifireads = Channel.fromPath(params.hifi)
                ontreads = Channel.fromPath(params.ont)
                input = genome.combine(hifireads)
                hifi = MINIMAP2(input, "hifi").mappings
            if (params.ont != false) {
                input = genome.combine(ontreads)
                ont = MAP_ONT(input, "ont").mappings
                alignments = hifi.concat(ont)}
            else{
                alignments = hifi
            }
        
            } else {
            alignments = Channel.fromPath(params.bam)}
        vcf = FREEBAYES_SINGLE(genome, hifi, 80, "").vcf
    } else {
        vcf = Channel.fromPath(params.vcf)
        hifibam = Channel.fromPath(params.hifibam)
        ontbam = Channel.fromPath(params.ontbam)
        alignments = hifibam.concat(ontbam).collect()
    }
    if (params.filter == true) {
        vcf2 = REFERENCE_VCF(vcf, [8,50]).vcf}
    else {
        vcf2 = vcf}
    haplotypes = PHASE_REFERENCE(genome, vcf2, alignments)
}
