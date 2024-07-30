include { BWA; BAMQC; FREEBAYES_JOINT; MARKDUP;FREEBAYES_GENOTYPE} from './modules/minimap.nf'  

include { JOIN_VCF; FILTER as DECOMPOSE } from './Popgen.nf'

include { PHASING } from './modules/phasing.nf'
params.winsize=10000

process FILTER{
    label 'low_resources'
    publishDir "${params.publishDir}/population/variant_calling/01_filtering", mode: 'copy'   
    label 'vcflib'
    maxForks 10

    input:
    path vcf
    tuple val(mindep), val(maxdep)

    output:
    path "*_filtered.vcf.gz", emit: vcf
    tuple path(".command.sh"), path(".command.log")

    script:
    """
    vcffilter -g 'DP > $mindep & DP < $maxdep' $vcf > ${vcf.simpleName}_filtered.vcf
    bgzip ${vcf.simpleName}_filtered.vcf
    """

}

process PLINK2_RAD {
    label 'low_resources'
    publishDir "${params.publishDir}/population/popgen/00_input", mode: 'copy'   

    input:
    path vcf
 
    output:
    path ".*command.*"
    tuple path ("*eigenvec"), path("*eigenval")
    path "${vcf.simpleName}_plink*"
    path "${vcf.simpleName}_plink*.vcf", emit: vcf

    
    script:

    """
    LC_NUMERIC="C" seq -s " " 0.05 0.05 1.0 > mybins.txt
    plink2 --make-pgen --vcf $vcf --allow-extra-chr --maf 0.05 --geno $params.missing_loci --export vcf-4.2 --out ${vcf.simpleName}  --set-all-var-ids '@:#\$r\$a' --new-id-max-allele-len 1000 truncate --threads $task.cpus --vcf-half-call missing
    plink2 --pfile ${vcf.simpleName} --allow-extra-chr --freq alt1bins-file='mybins.txt' --out ${vcf.simpleName}_afs
    plink2 --pfile ${vcf.simpleName} --allow-extra-chr --snps-only --read-freq ${vcf.simpleName}_afs.afreq --export vcf-4.2 --pca --out ${vcf.simpleName}_plink --maf 0.05    
    mv .command.sh .${vcf.simpleName}.command.sh
    mv .command.log .${vcf.simpleName}.command.log
    """
// missingness // hw 
}

process MERGEVCF{
    label 'vcflib'
    input:
    
    tuple path(vcf1), path(vcf2)
    output:
    path ("merged.vcf.gz"), emit: vcf

    script:
    """
tabix $vcf1
tabix $vcf2
vcf-merge $vcf1 $vcf2 | bgzip -c > merged.vcf.gz"""

}
process VCFtoBED{
    label 'vcflib'
    input:
    path vcf
    output:
    path "*bed", emit: bed
    script:
    """vt view $vcf | vcf2bed.py - | sort -k1,1 -k2,2n - > ${vcf.simpleName}.bed"""
}
params.aligned=false
params.missing_loci = 0.2
params.coverage = 3000
params.vcf =false
params.phasing = false 
params.merge_vcf = false

workflow {
// prepare files
    ASSEMBLY = Channel.fromPath(params.genome)    
    READS = channel.fromFilePairs(params.reads, checkIfExists:true)
// map and process
    if (params.aligned == false){
        INPUT = READS.combine(ASSEMBLY)
        MAPPINGS = BWA(INPUT).mappings
        NODUP = MARKDUP(MAPPINGS).mappings
        BAMQC(NODUP)
    } else {
        NODUP = Channel.fromPath(params.aligned, checkIfExists:true)     
    }
       
    COLLECTED_MAPPINGS = NODUP.collect()

// prepare reference
    contigs_file = file(params.contigs)
    allcontigs = Channel.fromList(contigs_file.readLines())
    GENOME_CONTIG = ASSEMBLY.combine(allcontigs)

// genotype based on vcf
    if (params.vcf == false){   
        VCFS = FREEBAYES_JOINT(GENOME_CONTIG, COLLECTED_MAPPINGS, params.coverage, "--use-duplicate-reads").vcf.collect()
        VCF_JOINED = JOIN_VCF(VCFS).vcf
        FILTERED = DECOMPOSE(VCF_JOINED, ASSEMBLY, [2, 50])
    } else {
        VCF = Channel.fromPath(params.vcf)
        bed = VCFtoBED(VCF).bed
        VCFSOUT = FREEBAYES_GENOTYPE(bed,ASSEMBLY, COLLECTED_MAPPINGS).vcf
        FILTERED = FILTER(VCFSOUT, [2,20])
        if (params.merge_vcf != false){
            INPUTVCF = VCF.combine(FILTERED.vcf) 
            FILTERED = MERGEVCF(INPUTVCF)
        }
    }    
    FINALVCF = PLINK2_RAD(FILTERED.vcf).vcf
    if (params.phasing != false){
        PHASED = PHASING(FINALVCF)
    }
}