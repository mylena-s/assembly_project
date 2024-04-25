// PUBLISHDIR????
process MINIMAP{
    publishDir "${params.publishDir}/population/alignment", mode: 'copy'   
    label 'medium_resources'
    label 'minimap'
    maxForks 10

    input:
    tuple path(genome), path(reads)

    output:
    path '*.bam', emit: mappings
    path ".command*"

    script:
    def mem = "${task.memory}".replaceAll("\\sGB","G")
    def name = "${reads.baseName}".replaceAll(".fastq","")
    """minimap2 -ax map-$params.dtype $genome $reads -o ${name}.sam
    samtools view -b ${name}.sam -o ${name}.bam 
    rm *sam
    bamaddrg -b ${name}.bam -s ${name} > ${name}_RG.bam 
    samtools sort -@${task.cpus} -O BAM -o ${name}_sorted.bam ${name}_RG.bam   
    rm ${name}.bam ${name}_RG.bam
    """
}

process MINIMAP2{
    label 'medium_resources'
    label 'minimap'
    maxForks 10
    publishDir "${params.publishDir}/alignment", mode: 'copy'   


    input:
    tuple path(genome), path(reads)
    val (type)

    output:
    path '*.bam', emit: mappings
    path '.*command.*'

    script:
    def mem = "${task.memory}".replaceAll("\\sGB","G")
    def name = "${reads.baseName}".replaceAll(".fastq","")

    """minimap2 -ax map-$type $genome $reads -o ${name}.sam
    samtools view -b ${name}.sam -o ${name}.bam 
    rm *sam
    samtools sort -@${task.cpus} -O BAM -o ${genome.baseName}_${name}_sorted.bam ${name}.bam   
    rm ${name}.bam
    mv .command.sh .${genome.baseName}_${name}.command.sh
    mv .command.log .${genome.baseName}_${name}.command.log
    
    """
}

process MARKDUP{
    label 'medium_resources'
    publishDir "${params.publishDir}/population/alignment", mode: 'copy'   

    input: 
    path bam

    output:
    path "*sorted_dup.bam", emit: mappings
    path ".command*"

    script:
    """sambamba markdup $bam ${bam.baseName}.sorted_dup.bam"""
}

process BAMQC {
    publishDir "${params.publishDir}/population/alignment/QC", mode: 'copy'    
    label 'low_resources'
    label 'qualimap'

    input:
    path bam
    
    output:
    path "${bam.baseName}/", emit: report
    path ".command*"

    script:
    """qualimap bamqc -bam $bam -nt $task.cpus -outdir ${bam.baseName} -outformat HTML
    """
}
process FREEBAYES_JOINT{
    label 'low_resources'
    publishDir "${params.publishDir}/population/SNP_variant_calling", mode: 'copy'   
    label 'minimap'
    maxForks 10

    input:
    tuple path(genome), val(sequence_id)
    val(bams)

    output:
    path "*_jointcall.vcf", emit: vcf
    path "jointcall.command.*"

    script:
    def files = bams.join(' ')
    """ls $files > bam_list
    freebayes -f $genome -L bam_list --gvcf -g 100 --ploidy 2 -r $sequence_id > ${genome.baseName}_${sequence_id}_jointcall.vcf
    mv .command.sh jointcall.command.sh
    mv .command.log jointcall.command.log
    """
}   

