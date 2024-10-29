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
    label 'resource_intensive'
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
    samtools view -bq 20 ${name}.sam -o ${name}.bam 
    rm *sam
    samtools sort -@${task.cpus} -O BAM -o ${genome.baseName}_${name}_sorted.bam ${name}.bam   
    rm ${name}.bam
    mv .command.sh .${genome.baseName}_${name}.command.sh
    mv .command.log .${genome.baseName}_${name}.command.log
    """
}


process BWA{
    label 'medium_resources'
    label 'minimap'
    maxForks 10
    publishDir  "${params.publishDir}/population/alignment", mode: 'copy'  

    input:
    tuple val(name), path(reads), path(genome)
    
    output:
    path '*_sorted.bam', emit: mappings
    path '.*command.*'

    script:
    def mem = "${task.memory}".replaceAll("\\sGB","G")

    """
    bwa index $genome
    bwa mem -t $task.cpus -R '@RG\\tID:${name}\\tSM:${name}' $genome ${reads[0]} ${reads[1]} | samtools view -bq 20  - > alignment.raw.bam
    samtools sort -@${task.cpus} -O BAM -o ${name}_sorted.bam alignment.raw.bam
    mv .command.sh .${name}.command.sh
    mv .command.log .${name}.command.log"""
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
    """sambamba markdup $bam ${bam.simpleName}.sorted_dup.bam"""
}

process BAMQC {
    publishDir "${params.publishDir}/population/alignment/QC", mode: 'copy'    
    label 'low_resources'
    label 'qualimap'

    input:
    path bam
    
    output:
    path "${bam.baseName}/", emit: report
    path "${bam.baseName}/genome_results.txt", emit: txt
    path ".command*"

    script:
    def mem = "${task.memory}".replaceAll("\\sGB","G")
    """qualimap bamqc -bam $bam --java-mem-size=$mem -nt $task.cpus -outdir ${bam.baseName} -outformat HTML
    """
}
process FREEBAYES_JOINT{
    label 'low_resources'
    publishDir "${params.publishDir}/population/variant_calling/00_freebayes", mode: 'copy'   
    label 'minimap'
    maxForks 5

    input:
    tuple path(genome), val(sequence_id)
    val(bams)
    val(cov)
    val(optional)

    output:
    path "*_jointcall.vcf", emit: vcf
    path "jointcall.command.*"

    script:
    def files = bams.join(' ')
    """
    samtools index -M $files
    ls $files > bam_list
    freebayes -f $genome -L bam_list --gvcf -g $cov --ploidy 2 -r $sequence_id $optional> ${genome.baseName}_${sequence_id}_jointcall.vcf
    mv .command.sh jointcall.command.sh
    mv .command.log jointcall.command.log
    """
}   

process FREEBAYES_GENOTYPE{
    label 'low_resources'
    publishDir "${params.publishDir}/population/variant_calling/00_freebayes_genotype", mode: 'copy'   
    label 'minimap'
    maxForks 5

    input:
    each path(bed)
    path(genome)
    val(bams)

    output:
    path "*_jointcall.vcf", emit: vcf
    path "jointcall.command.*"

    script:
    def files = bams.join(' ')
    """ls $files > bam_list
    freebayes -f $genome -t $bed -L bam_list --gvcf --ploidy 2 --report-monomorphic > ${genome.baseName}_regenotyped_jointcall.vcf
    mv .command.sh jointcall.command.sh
    mv .command.log jointcall.command.log
    """
}


process FREEBAYES_SINGLE{
    label 'low_resources'
    publishDir "${params.publishDir}/phasing/variant_calling/00_freebayes", mode: 'copy'   
    label 'vcflib'
    maxForks 5

    input:
    path(genome)
    val(bam)
    val(cov)
    val(optional)

    output:
    path "*.vcf", emit: vcf
    path ".command.*"

    script:
    """
    samtools index -M $bam
    freebayes -f $genome -b $bam --gvcf -g $cov --ploidy 2 > ${genome.baseName}.vcf
    mv .command.sh jointcall.command.sh
    mv .command.log jointcall.command.log
    """
}   