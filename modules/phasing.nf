process FIXVCF{
    label 'minimap'
    input:
    path(vcf)
    output:
    path("fixed.vcf"), emit: vcf
    script:
    """bcftools +fixploidy $vcf > fixed.vcf"""
}

process PHASING {
    label 'resource_intensive'
    label 'twiss'
    publishDir "${params.publishDir}/population/phylogeny/00_phasing/" , mode: 'copy'

    input:
    path(vcf)
    output:
    path("${vcf.baseName}.geno.gz"), emit:geno
    path("phased.vcf.gz")
    path(".command.*")

    script:
    """
    set JAVA_OPTS = "-Xmx12g"
    beagle.r1399.jar gt=$vcf out=phased impute=true nthreads=$task.cpus window=10000 overlap=1000 gprobs=false
    parseVCF.py -i phased.vcf.gz --skipIndels --minQual 20 --gtf flag=DP min=1 | gzip > ${vcf.simpleName}.geno.gz"""
}

process PHASE_VCF{
    label 'whatshap'
    label "resource_intensive"
    publishDir "${params.publishDir}/population/phylogeny/00_phasing/", mode: 'copy'

    input:
    path (reference)
    path(vcf)
    path(bams)

    output:
    path("phased.vcf"), emit: phased_vcf
    path("*phasing_statistics.tsv")
    tuple path(".whatshap.command.sh"), path(".whatshap.command.log")

    script:
    """
    samtools faidx $reference
    samtools index -M $bams -@$task.cpus
    whatshap phase -o phased.vcf --reference=$reference $vcf $bams
    for sample in *_sorted.sorted_dup.bam; do f=\$(basename -s "_sorted.sorted_dup.bam" "\$sample" ); whatshap stats phased.vcf --tsv \$f.phasing_statistics.tsv --sample \$f; done
    mv .command.sh .whatshap.command.sh
    mv .command.log .whatshap.command.log
    """
//    whatshap phase -o phased.vcf --merge-reads --error-rate=0.03 --maximum-error-rate=0.05 --chromosome Ctuca_Chr1 --reference=$reference $vcf $bams

}
