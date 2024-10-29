include { BWA; BAMQC; FREEBAYES_JOINT; MARKDUP;FREEBAYES_GENOTYPE} from './modules/minimap.nf'  

include { JOIN_VCF; FILTER as DECOMPOSE; VCF_QC } from './Popgen.nf'

include { PHASING } from './modules/phasing.nf'
include { PHYLONETWORK; TOPOLOGY_WEIGHTS } from './Phylogeny.nf'

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
    tuple path ("*pgen"), path("*psam"), path("*pvar")
    path "${vcf.simpleName}_pca.vcf", emit: vcf

    
    script:

    """
    plink2 --make-pgen --vcf $vcf --allow-extra-chr --geno 0.5 --out ${vcf.simpleName}_filter1 --set-all-var-ids '@:#\$r\$a' --new-id-max-allele-len 1000 truncate --threads $task.cpus --vcf-half-call missing
    plink2 --make-pgen --pfile ${vcf.simpleName}_filter1 --allow-extra-chr --snps-only --mind 0.7 --geno $params.missing_loci --freq --out ${vcf.simpleName}_filter2 --threads $task.cpus 
    plink2 --pfile ${vcf.simpleName}_filter2 --allow-extra-chr --read-freq ${vcf.simpleName}_filter2.afreq --export vcf-4.2 --pca --out ${vcf.simpleName}_pca    
    """
// missingness // hw 
}

process MERGEVCF{
    label 'minimap'
    input:
    
    tuple path(vcf1), path(vcf2)
    output:
    path ("merged.vcf.gz"), emit: vcf

    script:
    """
tabix $vcf1
tabix $vcf2
bcftools merge $vcf1 $vcf2 | bgzip -c > merged.vcf.gz"""

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
params.regenotype = true
params.phylo = false
params.variantcall = true


workflow VARIANTCALL {
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
        if (params.regenotype==true){
            bed = VCFtoBED(VCF).bed
            VCFSOUT = FREEBAYES_GENOTYPE(bed,ASSEMBLY, COLLECTED_MAPPINGS).vcf
            FILTERED = FILTER(VCFSOUT, [2,20])
            if (params.merge_vcf != false){
                INPUTVCF = VCF.combine(FILTERED.vcf) 
                FILTERED = MERGEVCF(INPUTVCF).vcf
            }
        } else {
            FILTERED = Channel.fromPath(params.regenotype)
            INPUTVCF = VCF.combine(FILTERED)
            FILTERED = MERGEVCF(INPUTVCF).vcf
        }
    }
    VCF_QC(FILTERED, ASSEMBLY)    
    FINALVCF = PLINK2_RAD(FILTERED).vcf

    if (params.phylo != false){
        PHASED_GENO = PHASING(FINALVCF).geno
        PHYLO(PHASED_GENO)}        
}

workflow PHYLO {
    take:
        PHASED_GENO
    main:
        DIST = PHYLONETWORK(PHASED_GENO)
        TWISST = TOPOLOGY_WEIGHTS(PHASED_GENO, ["3","5","10"])
}

workflow {
    if (params.variantcall == true) {
        VARIANTCALL()}
    if (params.phylo != false) {
        PHASED_GENO = Channel.fromPath(params.phased)
        PHYLO(PHASED_GENO)
    }

} 