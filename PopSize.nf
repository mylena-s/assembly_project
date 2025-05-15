process MSMC2{
    label 'low_resources'
    label 'minimap'
    publishDir "${params.publishDir}/population/size", mode: 'copy'   

    input:
    tuple path(genome), path(bam), path(repeats), path(mask)
    val(meandep)
    output:
    path ('*multihetsep.txt')
    path ('*out.vcf.gz')
    path ("*out_mask.bed.gz")
    path ("*out.msmc2*")
    script:
    """
    samtools index $bam
    samtools faidx $genome
    bcftools mpileup -q 30 -Q 20 -C 50 -f $genome $bam -Ou | bcftools call -c -V indels | bamCaller.py $meandep out_mask.bed.gz | gzip -c > out.vcf.gz
    generate_multihetsep.py --mask out_mask.bed.gz --mask $mask --negative_mask $repeats out.vcf.gz > multihetsep.txt
    msmc2_Linux -t 6 -p -o out.msmc2 multihetsep.txt
"""
    }

process MSMC2_DEEPVARIANT{
    label 'low_resources'
    label 'minimap'
    publishDir "${params.publishDir}/population/size", mode: 'copy'   

    input:
    tuple path(vcf), path(repeats), path(mask), val(chr)

    output:
    path ('*multihetsep.txt')
    path ("*out.msmc2*")

    script:
    """
    grep $chr $mask > ${chr}.${mask}
    tabix $vcf
    bcftools view $vcf -r $chr > ${chr}.${vcf}
    bgzip ${chr}.${vcf}
    generate_multihetsep.py --mask ${chr}.${mask} --negative_mask $repeats ${chr}.${vcf}.gz > ${chr}.multihetsep.txt
    msmc2_Linux -t 6 -p '15*1' -o ${chr}.out.msmc2 ${chr}.multihetsep.txt
"""
}



process CALLABLEMASK { 
    label 'minimap'
    publishDir "${params.publishDir}/population/size", mode: 'copy'   

    input:
    path(bam)
    val(mindep)
    val(maxdep)

    output:
    path ('*_mask.bed'), emit: callablemask
    script:
    """
    export MOSDEPTH_Q0=NO_COVERAGE
    export MOSDEPTH_Q1=CALLABLE
    export MOSDEPTH_Q2=REPEAT
    samtools index $bam
    mosdepth --quantize 0:$mindep:$maxdep: callable $bam -n
    gunzip *bed.gz
    grep -v NO_COVERAGE *bed | awk '{NF--; print}' OFS='\t' > ${bam.simpleName}_mask.bed
    """
}


params.wsize = 100000
params.maxDepth = 40
params.valQual = 20
params.meandep = 30
params.mindep = params.meandep /2
params.maxdep = (params.meandep /2) + params.meandep
params.vcf = false

workflow {
    ASSEMBLY = Channel.fromPath(params.genome)    
    BAM = Channel.fromPath(params.bam)
    REPEATS = Channel.fromPath(params.satellites)
    MASK = CALLABLEMASK(BAM,params.mindep,params.maxdep).callablemask
    if (params.vcf == false){
        INPUT = ASSEMBLY.combine(BAM).combine(REPEATS).combine(MASK)
        SIZE = MSMC2(INPUT, params.meandep)}
    if (params.vcf != false){
        VCF = Channel.fromPath(params.vcf)
        contigs_file = file(params.contigs)
        allcontigs = Channel.fromList(contigs_file.readLines())
        INPUT2 = VCF.combine(REPEATS).combine(MASK).combine(allcontigs)
        SIZ = MSMC2_DEEPVARIANT(INPUT2)
    }
}