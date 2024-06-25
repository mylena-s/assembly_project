
log.info """\
    P I P E L I N E     F O R    P O P G E N  
    ==========================================
   
    performs long read mapping with minimap
    sorts, marks duplicates, and calls variants    
    Current project dir. is: $projectDir"""
    .stripIndent()

process JOIN_VCF {
    label 'low_resources'
    label 'minimap'
    publishDir "${params.publishDir}/population/variant_calling/00_freebayes", mode: 'copy'   

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
    tabix -p vcf $vcf
    picard CreateSequenceDictionary -R $genome"""
}

process FILTER{
    label 'low_resources'
    publishDir "${params.publishDir}/population/variant_calling/01_filtering", mode: 'copy'   
    label 'vcflib'
    maxForks 10

    input:
    path vcf
    each path(reference)

    output:
    path "*_filtered_normalized.vcf.gz", emit: vcf
    tuple path(".command.sh"), path(".command.log")

    script:
    def name = "${vcf}".replaceAll(".vcf.gz","")
    """
    bgzip -d $vcf
    vt decompose_blocksub ${vcf.baseName} | vt normalize -r $reference - > intermediate.vcf
    vcffilter -f 'QUAL / AO > 10 & TYPE = snp' -g 'DP > 1 & DP < 20' intermediate.vcf > ${name}_filtered_normalized.vcf
    bgzip *_filtered_normalized.vcf
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


process PLINK2_FILTER {
    label 'low_resources'
    publishDir "${params.publishDir}/population/popgen/00_input", mode: 'copy'   

    input:
    path vcf

    output:
    path ".*command.*"
    path ("*vcf"), emit: vcf

    script:
    def name = "${vcf}".replaceAll(".vcf.gz","")
    """
    plink2 --make-pgen --vcf $vcf --allow-extra-chr --max-alleles 2 --mind 0.1 --export vcf-4.2 --vcf-half-call m --maf 0.05 --out $name --set-all-var-ids '@:#\$r\$a' --new-id-max-allele-len 1000 truncate
    mv .command.sh .${name}.command.sh
    mv .command.log .${name}.command.log
    """
    // no multiallelic, remove samples with > 10% missing, remove alelles with freq <5%
}

process PLINK2_FORMAT {
    label 'low_resources'
    publishDir "${params.publishDir}/population/popgen/00_input", mode: 'copy'   

    input:
    path vcf
    path phenotypes

    output:
    path ".*command.*"
    tuple path("*.pgen"), path("*.psam"), path("*.pvar"), emit: pfiles
    path ("*vcf"), emit: vcf

    script:
    def name = "${vcf}".replaceAll(".vcf.gz","")
    """
    plink2 --make-pgen --vcf $vcf --allow-extra-chr --max-alleles 2 --mind 0.1 --export vcf-4.2 --vcf-half-call m --maf 0.05 --out $name --set-all-var-ids '@:#\$r\$a' --new-id-max-allele-len 1000 truncate --pheno $phenotypes 
    mv .command.sh .${name}.command.sh
    mv .command.log .${name}.command.log
    """
    // no multiallelic, remove samples with > 10% missing, remove alelles with freq <5%
}

process PLINK2_WIDE {
    label 'low_resources'
    publishDir "${params.publishDir}/population/popgen/00_input", mode: 'copy'   

    input:
    tuple path(pgen), path(psam), path(pvar)
 
    output:
    path ".*command.*"
    tuple path("*.pgen"), path("*.psam"), path("*.pvar"), emit: pfiles
    path "*wide.vcf", emit: vcf

    script:
    def name = "${pgen.baseName}_wide"

    """
    plink2 --pfile ${pgen.baseName} --allow-extra-chr --maf 0.05 --geno 0.01 --export vcf-4.2 --out ${name} --threads $task.cpus --make-pgen
    mv .command.sh .${name}.command.sh
    mv .command.log .${name}.command.log
    """
// missingness // hw 
}

process PLINK2_PAIRWISE {
    label 'low_resources'
    publishDir "${params.publishDir}/population/popgen/00_input", mode: 'copy'   

    input:
    tuple path(pgen), path(psam), path(pvar)
    path population

    output:
    path ".*command.*"
    tuple path("*.pgen"), path("*.psam"), path("*.pvar"), emit: pfiles
    path "*_pairwise.vcf", emit: vcf

    script:
    def name = "${pgen.baseName}_pairwise"

    """
    plink2 --pfile ${pgen.baseName} --allow-extra-chr --maf 0.05 --geno 0.05 --threads $task.cpus --keep $population --hwe 1e-20 --keep-if morphology == 1 --write-snplist --out hwe_pass
    plink2 --pfile ${pgen.baseName} --extract hwe_pass.snplist --export vcf-4.2 --out ${name} --make-pgen
    mv .command.sh .${name}.command.sh
    mv .command.log .${name}.command.log
    """
// missingness // hw 
}


process AFS{
    label 'low_resources'
    publishDir "${params.publishDir}/population/popgen/02_afs", mode: 'copy'   

    input:
    tuple path(pgen), path(psam), path(pvar)
    each path(population)

    output:
    tuple path("*.pgen",includeInputs:true), path("*.psam",includeInputs:true), path("*.pvar",includeInputs:true), path("*.afreq"), path("*.bins"), emit: pfiles
    tuple path(".*command.sh"), path(".*command.log")
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
    publishDir "${params.publishDir}/population/popgen/02_afs", mode: 'copy'   
    label 'python'

    input:
    tuple path("*.pgen"), path("*.psam"), path("*.pvar"), path(afs), path("*.bins")

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
    publishDir "${params.publishDir}/population/popgen/03_fst", mode: 'copy'
    errorStrategy 'ignore'

    input:
    each path(vcf)
    tuple path (pop1), path(pop2)

    output:
    path "*.fst"
    path "*list"
    path ".*command.log"
    path ".*command.sh"
    
    script:
    def name = "${vcf.baseName}_${pop1.baseName}_${pop2.baseName}"
    """
    vcftools --gzvcf $vcf --weir-fst-pop $pop1 --weir-fst-pop $pop2 --fst-window-size $params.winsize --out ${name}_w$params.winsize
    awk 'NR==1 || \$5 > 0.4' *.fst > ${name}_highWeightedFst.list
    awk 'NR==1 || \$6 > 0.2' *.fst > ${name}_highMeanFst.list 
    mv .command.sh .${name}.command.sh
    mv .command.log .${name}.command.log
    """
}


process STATISTICS_WINDOWS{
    label 'low_resources'
    publishDir "${params.publishDir}/population/popgen/01_stats", mode: 'copy'
    errorStrategy 'ignore'

    input:
    path vcf
    each (pop)
    output:
    path "*.het"
    path "*.imiss"
    path "*.lmiss"
    path "*.pi"
    path ".*command.log"
    path ".*command.sh"
    path "*warnings.txt"
    
    script:
    def name = "${vcf.baseName}_${pop.baseName}"

    """
    vcftools --gzvcf $vcf --out ${name}_w$params.winsize --keep $pop --window-pi $params.winsize  2> warnings.txt
    vcftools --gzvcf $vcf --het --out ${name} --keep $pop 2>> ${name}_warnings.txt
    vcftools --gzvcf $vcf --missing-indv --out ${name} --keep $pop 2>> ${name}_warnings.txt
    vcftools --gzvcf $vcf --missing-site --out ${name} --keep $pop 2>> ${name}_warnings.txt
    mv .command.sh .${name}_command.sh
    mv .command.log .${name}_command.log
    """

}

process PLINK_PCA {
    publishDir "${params.publishDir}/population/popgen/04_pca", mode: 'copy'
    label 'low_resources'
    input:
    tuple path(pgen), path(psam), path(pvar), path(freq), path("*.bins")

    output:
    tuple path("*.eigenvec"), path("*.eigenval"),path(".*command.log"), emit:results
    tuple path(".*command.sh"), path(".*command.log")

    script:
    def name = "${pgen.baseName}"
    """
    plink2 --pfile $name --allow-extra-chr --snps-only --read-freq $freq --pca --out $name --maf 0.05
    mv .command.sh .${name}.command.sh
    mv .command.log .${name}.command.log
    """
}

process GWAS{
    label 'low_resources'
    publishDir "${params.publishDir}/population/popgen/05_gwas", mode: 'copy'

    input:
    tuple path(pgen), path(psam), path(pvar)
    output:
    path("gwas_results*")
    tuple path (".command.sh"), path(".command.log")

    script:
    """
    plink2 --pfile ${pgen.baseName} --pheno-name $params.phenotype_name --glm no-firth allow-no-covars --out gwas_results --allow-extra-chr --threads $task.cpus  
    """
//--covar ${covariateFile} --covar-col-nums ${covariateCols}  hide-covar https://cloufield.github.io/GWASTutorial/06_Association_tests/
} 

process PLOT_PCA {
    label 'python'
    publishDir "${params.publishDir}/population/popgen/04_pca", mode: 'copy'

    input:
    tuple path(eigenvec), path(eigenval), path(prevlog)
    path(sample_info)

    output:
    path("*PCA.svg")

    script:
    """
#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re

eigenval = pd.read_csv('${eigenval}', sep='\\s+', header=None)
eigenval["percent"] = eigenval[0] / sum(eigenval[0]) * 100
PC1 = "PC1 ("+str(eigenval.percent[0].round(1))+"%)"
PC2 = "PC2 ("+str(eigenval.percent[1].round(1))+"%)"
PC3 = "PC3 ("+str(eigenval.percent[2].round(1))+"%)"
pattern = re.compile("[0-9]+ variants remaining")

file = open('${prevlog}')
for i, line in enumerate(file):
    for match in re.finditer(pattern, line):
        n_SNPS = round(int(match[0].split(" ")[0])/1000000,1)

file.close()
fig, axis = plt.subplots(nrows=2, ncols=2, figsize=(10,10),)
screplot = axis[0,0].bar(x=eigenval.index+1 , height = eigenval.percent,
                         alpha=0.5)
axis[0,0].set_xlabel("Principal component", fontsize = 12)
axis[0,0].set_ylabel("Percentaje of explained variance (%)", fontsize = 12)
eigenvec = pd.read_csv('${eigenvec}', sep='\\t')
population = pd.read_csv("${sample_info}", sep="\\t")
eigenvec = pd.merge(eigenvec, population, left_on="#IID", right_on="#IID",
                    how="left")
PC1_2= sns.scatterplot(x=eigenvec.PC1, y=eigenvec.PC2, hue=eigenvec.morphology,
                       ax= axis[0,1],style=eigenvec.population, alpha=0.5,
                       legend='full', s=400)
legend=axis[0,1].legend(fontsize="large", bbox_to_anchor=(1.05, 1),
                        loc='upper left', title="References\\n ("+ str(n_SNPS)
                        + "M SNPs)")
axis[0,1].set_xlabel(PC1, fontsize = 12)
axis[0,1].set_ylabel(PC2, fontsize = 12)
plt.setp(legend.get_title(),fontsize='large')
PC1_3= sns.scatterplot(x=eigenvec.PC1, y=eigenvec.PC3, hue=eigenvec.morphology,
                       ax= axis[1,0],style=eigenvec.population, alpha=0.5,
                       s=400, legend=False)
axis[1,0].set_xlabel(PC1, fontsize = 12)
axis[1,0].set_ylabel(PC3, fontsize = 12)
PC3_3=sns.scatterplot(x=eigenvec.PC2, y=eigenvec.PC3, hue=eigenvec.morphology,
                      ax= axis[1,1],style=eigenvec.population, alpha=0.5,
                      s=400, legend=False)
axis[1,1].set_xlabel(PC2, fontsize = 12)
axis[1,1].set_ylabel(PC3, fontsize = 12)
plt.savefig("${eigenval.baseName}_PCA.svg", dpi=300)
"""
}

params.winsize=10000
reference = "$projectDir/lib/NC_030272.1"
params.assembled_mito = false
params.vcfs = false
params.phenotype = "$projectDir/lib/dummy.file"
params.phenotype_name = "morphology"
params.dtype = "hifi"
params.mapped = false
params.nodup = false

include { MINIMAP; BAMQC; FREEBAYES_JOINT; MARKDUP } from './modules/minimap.nf'  
include { PHASE_VCF; FIXVCF } from './modules/phasing.nf'  


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
        MAPPINGS = Channel.fromPath(params.mapped) 
    }
    if (params.nodup == false) {
        MAPPINGS_NODUP = MARKDUP(MAPPINGS).mappings
        COLLECTED_MAPPINGS = MAPPINGS_NODUP.collect()
    } else {
        MAPPINGS_NODUP = Channel.fromPath(params.nodup, checkIfExists:true) 
        COLLECTED_MAPPINGS = MAPPINGS_NODUP.collect()}
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
        ASSEMBLY
    main:
        INPUT = PREPARE_FILESQC(ASSEMBLY, VCF) 
        STATISTICS(INPUT)
}


workflow PAIRWISE{
    take:
    PFILES
    main:
        POPULATION1 = Channel.fromPath(params.population1, checkIfExists: true)
        POPULATION2 = Channel.fromPath(params.population2, checkIfExists: true)
        COMPARISONS = POPULATION1.combine(POPULATION2)
        PLINKOUTPUT = PLINK2_PAIRWISE(PFILES, COMPARISONS)
        PFILES = PLINKOUTPUT.pfiles 
        PLINK_VCF = PLINKOUTPUT.vcf  
        FST_WINDOWS(PLINK_VCF, COMPARISONS)
        GWAS(PFILES)
        //FREQSPECTRUM = AFS(PFILES, POPULATION) 
        //PLOTAFS(FREQSPECTRUM.pfiles)

}

workflow WIDE{
    take:
    PFILES
    PHENOTYPES
    main:
        POPULATION = Channel.fromPath(params.population, checkIfExists:true)
        PLINKOUTPUT = PLINK2_WIDE(PFILES)
        PFILES = PLINKOUTPUT.pfiles      
        PVCF=PLINKOUTPUT.vcf
        STATISTICS_WINDOWS(PVCF, POPULATION)
        FREQSPECTRUM = AFS(PFILES, POPULATION)    
        PCA = PLINK_PCA(FREQSPECTRUM.pfiles)
        PLOT_PCA(PCA.results, PHENOTYPES) //esto puede ser mejorado..d. en PCA.results ya podría estar la info de fenotipos... ya que esa info esta en psam
}


params.filter = true
params.vcf = false
params.QC = true
params.population = false
params.pca = false
params.pairwise = false
params.wide = false 


// main workflow following

workflow {
    ASSEMBLY = Channel.fromPath(params.genome)    

    if (params.vcf == false) {
        READS = Channel.fromPath(params.reads, checkIfExists:true)
        VCF = NUCLEAR(READS)
        // MITO(READS)
    }
    else {
        VCF = Channel.fromPath(params.vcf, checkIfExists: true)
    }

    if (params.filter == true) {
        filtered = FILTER(VCF, ASSEMBLY).vcf
    } else {
        filtered = VCF
    }

    if (params.QC == true){
        VCFs = VCF.concat(filtered)
        VCF_QC(VCFs, ASSEMBLY)
    } else {
        println("QC was skipped") 
    }

    VCFs = filtered
    PHENOTYPES = Channel.fromPath(params.phenotype, checkIfExists: true)
    PLINKVCF = PLINK2_FILTER(VCFs).vcf
    // phasing
    BAMS = Channel.fromPath(params.bams, checkIfExists:true).collect()
    VCF = FIXVCF(PLINKVCF).vcf
    PHASED_GENO = PHASE_VCF(ASSEMBLY, VCF, BAMS).phased_vcf

    //produce pfiles
    PLINKOUTPUT = PLINK2_FORMAT(PHASED_GENO, PHENOTYPES)
    PFILES = PLINKOUTPUT.pfiles

    if (params.pairwise != false) {
        PAIRWISE(PFILES)    
    } 
    if (params.wide != false) {
        WIDE(PFILES, PHENOTYPES)
    }
}
