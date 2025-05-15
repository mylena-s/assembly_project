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
    bcftools concat $files > ${sample}_all_joined.vcf
    bgzip ${sample}_all_joined.vcf
    mv .command.sh .concat.command.sh
    mv .command.log .concat.command.log
    """
}

params.technology = "PACBIO"
workflow {
//// remove from here
ASSEMBLY = Channel.fromPath(params.genome)    
MAPPINGS = Channel.fromPath(params.mapped)

//// to here
contigs_file = file(params.contigs)
allcontigs = Channel.fromList(contigs_file.readLines())
GENOME_CONTIG = ASSEMBLY.combine(allcontigs)
INPUT = GENOME_CONTIG.combine(MAPPINGS)
VCFS = DEEPVARIANT(INPUT, params.technology)
GROUPS = VCFS.genomicvcf.groupTuple(by:2)
SAMPLEVCF = JOIN_GVCF(GROUPS)
INPUT2 = SAMPLEVCF.vcf.collect()
INPUT2.count().subscribe { length -> if (length > 1) {
    GLXNEXUS(INPUT2)}}
}
