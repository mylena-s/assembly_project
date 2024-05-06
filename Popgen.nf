
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
    publishDir "${params.publishDir}/population/SNP_variant_calling/00_freebayes", mode: 'copy'   

    input:
    path(vcfs)

    output:
    path "all_jointcall.vcf.gz", emit: vcf
    path ".concat.command.*"

    script:
    def files = vcfs.join(' ')
    """
    bcftools concat $files > all_jointcall.vcf
    bgzip all_jointcall.vcf
    mv .command.sh .concat.command.sh
    mv .command.log .concat.command.log
    """
}

process STATISTICS{
    label 'low_resources'
    publishDir "${params.publishDir}/population/variant_calling/02_statistics", mode: 'copy'   
    label 'discvrseq'
    maxForks 10

    input:
    tuple path(genome), path(fai), path(dict), path(vfcfile), path(tbi) 

    output:
    path "*variantQC.html"
    
    script:
    """
    java -jar $projectDir/software/DISCVRSeq-1.3.72.jar VariantQC -R $genome -V $vfcfile -O ${vfcfile.baseName}.variantQC.html
    """
}

process PREPARE_FILESQC{
    label 'vcflib'
    label 'low_resources'

    input:
    path genome
    each path(vcf)

    output:
    tuple path("*fasta",includeInputs:true), path("*fasta.fai"), path("*dict"), path("*vcf.gz",includeInputs:true), path("*tbi") 

    script:
    """samtools faidx $genome
    tabix $vcf
    picard CreateSequenceDictionary -R $genome"""
}

process FILTER{
    label 'low_resources'
    publishDir "${params.publishDir}/population/variant_calling/01_filtering", mode: 'copy'   
    label 'vcflib'
    maxForks 10

    input:
    path vcf

    output:
    path "*_filtered.vcf.gz", emit: vcf

    script:
    def name = "${vcf}".replaceAll(".vcf.gz","")
    """
    tabix $vcf
    vcffilter -f 'QUAL > 20' -g 'DP > 1 & DP < 30' $vcf > ${name}_filtered.vcf
    bgzip *filtered.vcf
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
    path "${name}.final_mitogenome.fasta", emit: genome
    path ".*command.*"


    script:
    def name = "${read.baseName}".replaceAll(".fastq","")
    """
    mitohifi.py -r $read -f ${reference}.fasta -g ${reference}.gb -t $task.cpus -a animal 
    sed 's/>/>$name /g' final_mitogenome.fasta > ${name}.final_mitogenome.fasta
    mv final_mitogenome.gb ${name}.final_mitogenome.gb
    mv .command.log .${name}.command.log
    mv .command.sh .${name}.command.sh
    """
}


process PLINK2_FORMAT {
    label 'low_resources'
    publishDir "${params.publishDir}/population/popgen/00_input", mode: 'copy'   

    input:
    path vcf
    output:
    path ".*command.*"
    tuple path("*.pgen"), path("*.psam"), path("*.pvar"), emit: pfiles

    script:
    def name = "${vcf}".replaceAll(".vcf.gz","")
    """
    plink2 --make-pgen --vcf $vcf --allow-extra-chr --vcf-half-call m --out $name
    mv .command.sh .${name}.command.sh
    mv .command.log .${name}.command.log
    """

}

process AFS{
    label 'low_resources'
    publishDir "${params.publishDir}/population/popgen/01_afs", mode: 'copy'   

    input:
    tuple path(pgen), path(psam), path(pvar)
    each path(population)

    output:
    path "*.afreq", emit: freq
    path "*.bins", emit: bins

    script:
    def name = "${pgen.baseName}_${population}".replaceAll(".pop","")
    """
    LC_NUMERIC="C" seq -s " " 0.05 0.05 1.0 > mybins.txt
    plink2 --pfile ${pgen.baseName} --allow-extra-chr --keep $population --freq alt1bins-file='mybins.txt' --out ${name}_afs
    mv .command.sh .${name}.command.sh
    mv .command.log .${name}.command.log
    """
}

process PLOTAFS{
    label 'low_resources'
    publishDir "${params.publishDir}/population/popgen/01_afs", mode: 'copy'   
    label 'python'

    input:
    path afs
    output:
    path '*png'

    script:
    """#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
df = pd.read_csv('${afs}', sep="\t", skiprows=1, header=None)
plt.bar(x=df[0], height=df[1])
plt.xlabel("Alternate allele frequency", fontsize = 10)
plt.ylabel("Number of loci", fontsize = 10)
plt.savefig("${afs.baseName}.png")"""
}


process FST_WINDOWS{
    label 'low_resources'
    publishDir "${params.publishDir}/population/popgen/02_fst", mode: 'copy'   

    input:
    path vcf
    tuple path (pop1), path(pop2)

    output:
    path "*.fst"
    path "*list"
    
    script:
    """
    vcftools --gzvcf $vcf --weir-fst-pop $pop1 --weir-fst-pop $pop2 --fst-window-size $params.winsize
    --out ${pop1.baseName}_${pop2.baseName}_w$params.winsize
    awk 'NR==1 || \$5 > 0.4' *.fst > highWeightedFst.list
    awk 'NR==1 || \$6 > 0.2' *.fst > highMeanFst.list 
    mv .command.sh .${pop1.baseName}_${pop2.baseName}.command.sh
    mv .command.log .${pop1.baseName}_${pop2.baseName}.command.log
    """
}
params.winsize=10000


reference = "$projectDir/lib/NC_030272.1"
params.assembled_mito = false
params.vcfs = false

workflow MITO {
    take:
    READS
    main:
    if (params.assembled_mito == false){
        MITOCONDRIAS = MITOHIFI(READS,reference).genome.collect()}
    else {
        MITOCONDRIAS = Channel.fromPath(params.assembled_mito).collect()}
//    alignments = ALIGN_MITO(MITOCONDRIAS)
     
    }
    
// MITO workflow is not working. Problems with multiple genome alignment
// for pop.genomics better to do variant calling

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
    BAMQC(MAPPINGS_NODUP)
    contigs_file = file(params.contigs)
    allcontigs = Channel.fromList(contigs_file.readLines())
    GENOME_CONTIG = ASSEMBLY.combine(allcontigs)
    if (params.vcfs == false) {
        VCFS = FREEBAYES_JOINT(GENOME_CONTIG, COLLECTED_MAPPINGS).vcf.collect()
    } else {
        VCFS = Channel.fromPath(params.vcfs, checkIfExists: true).collect()}   
    ALL_VCF = JOIN_VCF(VCFS).vcf
    emit:
    ALL_VCF
}

workflow VCF_QC{
    take:
        VCF
    main:
        ASSEMBLY = Channel.fromPath(params.genome, checkIfExists: true)

    INPUT = PREPARE_FILESQC(ASSEMBLY, VCF) 
    STATISTICS(INPUT)
}

params.filter = true
params.vcf = false
params.QC = true
params.population = false

workflow {
    if (params.vcf == false) {
        READS = Channel.fromPath(params.reads, checkIfExists:true)
        VCF = NUCLEAR(READS)
        // MITO(READS)
    }
    else {
        VCF = Channel.fromPath(params.vcf, checkIfExists: true)
    }

    if (params.filter == true) {
        filtered = FILTER(VCF).vcf
    } else {
        filtered = channel.empty()
    }

    VCFs = VCF.concat(filtered)

    if (params.QC == true){
        VCF_QC(VCFs)
    } else {
        println("QC was skipped") 
    }

    if (params.population != false){
        POPULATION = Channel.fromPath(params.population, checkIfExists: true)
        PLINKOUTPUT = PLINK2_FORMAT(VCFs).pfiles
        FREQSPECTRUM = AFS(PLINKOUTPUT, POPULATION)}
        PLOTAFS(FREQSPECTRUM.bins)
        POPULATION2 = POPULATION
        COMPARISONS = POPULATION.combine(POPULATION2)
        FST_WINDOWS(VCFs, COMPARISONS)
    }
    
