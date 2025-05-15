params.vcf = false
params.wsize = 100000
params.meandep = 30
params.minDepth = params.meandep /2
params.maxDepth = (params.meandep /2) + params.meandep
params.valQual = 20
params.hom = false
params.map = false
params.call = false
params.technology = "PACBIO"
params.dtype = "hifi"
params.variant_analysis = false
params.popsize = false
params.bam = false 

process HETEROZYGOSIS {
    label 'low_resources'
    label 'minimap'
    publishDir "${params.publishDir}/population/heterozygosis", mode: 'copy'   

    input:
    tuple path(vcf), path(bams),path(windows)
    val(depth)
    val(vq)

    output:
    path("*bed"), emit: density
    path("*.het.vcf"), emit: hetsites
    path ".command.*"

    script:
    """
    bcftools filter -i 'QUAL > $vq' $vcf | bcftools filter -i 'FMT/DP < $depth' | bcftools view -v snps | bcftools view -g het > ${vcf.simpleName}_filtered.het.vcf 
    bedtools coverage -a windows.bed -b ${vcf.simpleName}_filtered.het.vcf > ${vcf.simpleName}_filtered.het.bed
    """

}



process HOMOZYGOSIS {
    label 'low_resources'
    label 'minimap'
    publishDir "${params.publishDir}/population/homozygosity", mode: 'copy'   

    input:
    tuple path(vcf), path(bams), path(windows)
    val(depth)
    val(mindepth)
    val(vq)

    output:
    path("*bed"), emit: density
    path("*.hom.vcf"), emit: homsites
    path ".command.*"

    script:
    """
    bcftools filter -i 'QUAL > $vq' $vcf | bcftools filter -i 'FMT/DP < $depth' | bcftools filter -i 'FMT/DP > $mindepth'| bcftools view -v snps | bcftools view -g hom > ${vcf.simpleName}_filtered.hom.vcf 
    bedtools coverage -a windows.bed -b ${vcf.simpleName}_filtered.hom.vcf > ${vcf.simpleName}_filtered.hom.bed"""
}
process MAKE_WINDOWS {
    label 'low_resources'
    label 'minimap'
    input:
    path genome
    val wsize

    output:
    path "windows.bed", emit: windows

    script:
    """samtools faidx $genome
    bedtools makewindows -g $genome'.fai' -w $wsize > windows.bed
    """
}


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
    msmc2_Linux -t 6 -o out.msmc2 multihetsep.txt
"""
    }

process MSMC2_DEEPVARIANT{
    label 'low_resources'
    label 'minimap'
    publishDir "${params.publishDir}/population/size", mode: 'copy'   

    input:
    tuple path(vcf), path(bam), path(repeats), path(mask), val(chr)

    output:
    path ('*multihetsep.txt')
    path ("*out.msmc2*")
    path ("*vcf.gz")

    script:
    """
    grep ${chr} ${mask} > ${chr}.${mask}
    tabix $vcf
    bcftools filter -i 'QUAL > $params.valQual' $vcf | bcftools filter -i 'FMT/DP < $params.maxDepth' | bcftools filter -i 'FMT/DP > $params.minDepth'| bcftools view -v snps | bcftools view -i 'strlen(ALT[0])==1 & strlen(ALT[1])<2' > ${vcf.simpleName}_filtered.vcf 
    bgzip ${vcf.simpleName}_filtered.vcf
    tabix ${vcf.simpleName}_filtered.vcf.gz
    bcftools view ${vcf.simpleName}_filtered.vcf.gz -r $chr > ${chr}.${vcf}
    bgzip ${chr}.${vcf}
    generate_multihetsep.py --mask ${chr}.${mask} --negative_mask $repeats ${chr}.${vcf}.gz > ${chr}.multihetsep.txt
    msmc2_Linux -t $task.cpus -p '15*1' -o ${chr}.out.msmc2 ${chr}.multihetsep.txt

"""
}

process CALLABLEMASK { 
    label 'minimap'
    publishDir "${params.publishDir}/population/size", mode: 'copy'   

    input:
    tuple path(vcf), path(bam)
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

process DEEPVARIANT{
publishDir "${params.publishDir}/population/variant_calling/deepvariant/", mode: 'copy'
maxForks 2

label 'deepvariant'
label 'resource_intensive'

input:
tuple path(genome), val(contig), path(mapping)
val(model)

output:
tuple path("*.gvcf"), val(contig), val("${mapping.simpleName}"), emit: genomicvcf
path("*.vcf"), emit: vcf
path("*.visual_report.html"), emit: report

script:
"""
samtools faidx $genome
samtools index $mapping
mkdir tmp
mkdir intermediate
/opt/deepvariant/bin/run_deepvariant --model_type $model --ref $genome --reads $mapping --output_vcf ${mapping.simpleName}_${contig}.vcf --output_gvcf ${mapping.simpleName}_${contig}.gvcf --regions $contig --num_shards $task.cpus --intermediate_results_dir intermediate
"""
}

//process SURVIVOR?{}
process GLXNEXUS {
    label 'minimap'
    label 'medium_resources'

    input:
    path (gvcfs)

    output:
    path("allmerged.bcf"), emit: bcf
    path("allmerged.vcf.gz")

    tuple path(".command.log"), path(".command.sh")

    script:
    """glnexus_cli --config DeepVariant *vcf.gz --threads $task.cpus > allmerged.bcf
    bcftools view allmerged.bcf | bgzip -@ $task.cpus -c > allmerged.vcf.gz
    """
} /// para q funcione este comando, cada gvcf tiene q tener el nombre de la muestra en lugar del nombre "default" para la columna genotipo
 // No puede hacer un único merging para todos los cromosomas si estos no estan mergeados
process JOIN_GVCF {
    label 'low_resources'
    label 'minimap'
    publishDir "${params.publishDir}/population/variant_calling/", mode: 'copy'   

    input:
    tuple path(vcfs), val(contigs), val(sample) 

    output:
    path("*all_joined.vcf.gz"), emit: vcf
    path ".concat.command.*"

    script:
    def files = vcfs.join(' ')
    """
    bcftools concat $files > ${sample}_all_joined.tmp.vcf
    sed 's/default/${sample}/g' ${sample}_all_joined.tmp.vcf > ${sample}_all_joined.vcf
    bgzip ${sample}_all_joined.vcf
    mv .command.sh .concat.command.sh
    mv .command.log .concat.command.log
    """
}

workflow CALL{
    take:
        ASSEMBLY
    
    main:
        if (params.bam == false) {
            READS = Channel.fromPath(params.reads)
            INPUT = ASSEMBLY.combine(READS)
            MAPPINGS = MINIMAP2(INPUT, params.dtype).mappings
        } else {
            MAPPINGS = Channel.fromPath(params.bam)
        }
        contigs_file = file(params.contigs)
        allcontigs = Channel.fromList(contigs_file.readLines())
        GENOME_CONTIG = ASSEMBLY.combine(allcontigs)
        INPUT2 = GENOME_CONTIG.combine(MAPPINGS)
        VCFS = DEEPVARIANT(INPUT2, params.technology)
        GROUPS = VCFS.genomicvcf.groupTuple(by:2)
        SAMPLEVCF = JOIN_GVCF(GROUPS)
        JOINED_VCFS = SAMPLEVCF.vcf.collect()
        VCF = JOINED_VCFS
            .map { it.size() }
            .ifElse(
                { count -> count > 1 },
                { GLXNEXUS(JOINED_VCFS).vcf },
                { JOINED_VCFS }
            )
    OUTPUT = VCF.combine(MAPPINGS)
    emit:
        OUTPUT
}


include { MINIMAP2} from './modules/minimap.nf'  

workflow VARIANT_ANALYSIS{
    take:
        VCF_BAMS
        ASSEMBLY
    main:
        WINDOWS = MAKE_WINDOWS(ASSEMBLY, params.wsize).windows
        HINPUT = VCF_BAMS.combine(WINDOWS)
        HETDENSITY = HETEROZYGOSIS(HINPUT, params.maxDepth, params.valQual)
        if (params.hom != false){
            HOMDENSITY = HOMOZYGOSIS(HINPUT, params.maxDepth,params.minDepth, params.valQual)
        }
}

workflow POPSIZE {
    take:
        VCF_BAM
    main:
        REPEATS = Channel.fromPath(params.satellites)
        MASK = CALLABLEMASK(VCF_BAM,params.minDepth,params.maxDepth).callablemask
        contigs_file = file(params.contigs)
        allcontigs = Channel.fromList(contigs_file.readLines())
        INPUT2 = VCF_BAM.combine(REPEATS).combine(MASK).combine(allcontigs)
        SIZ = MSMC2_DEEPVARIANT(INPUT2)
    }

workflow {
    ASSEMBLY = Channel.fromPath(params.genome)
    if (params.call != false){
        INPUT = CALL(ASSEMBLY)        
    } else {
        VCF = Channel.fromPath(params.vcf) /// es importante q este vcf sea vcf y no gvcf, de la manera q está escrito el pipeline aquí va a dar error, por lo tanto hay q correr en dos etapas: calling y luego todo el resto
        BAM = Channel.fromPath(params.bam)
        INPUT = VCF.combine(BAM)} 
    if (params.variant_analysis != false) {
        VARIANT_ANALYSIS(INPUT, ASSEMBLY)
    }    
    if (params.popsize != false ) {
        POPSIZE(INPUT)
    }

}