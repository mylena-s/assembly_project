process ESTIMATE_TREES{
    label 'low_resources'
    label 'twiss'
    publishDir "${params.publishDir}/population/phylogeny/01_topologyWeights/" , mode: 'copy'

    input:
    path(phased_genotypes)
    each size 

    output:
    path("phased.geno.phyml_bionj.w*.trees.gz"), emit: trees
    path("phased.geno.phyml_bionj.w*.tsv")
    tuple path(".phyml.command.log"),path(".phyml.command.sh")

    script:
    """
    parseVCF.py -i $phased_genotypes --skipIndels > phased.geno                                                                    
    phyml_sliding_windows.py -T $task.cpus -g phased.geno --prefix phased.geno.phyml_bionj.w${size} -w $size --windType sites --model GTR --optimise n
    mv .command.sh .phyml.command.sh
    mv .command.log .phyml.command.log"""
}

process TWISST {
    label 'medium_resources'
    publishDir "${params.publishDir}/population/phylogeny/01_topologyWeights/" , mode: 'copy'

    label 'twiss'

    input:
    path(trees)
    each path(groups)

    output:
    path("*.weights.csv.gz")
    tuple path(".twisst.command.sh"), path(".twisst.command.log")

    
    script:
    """
    twisst.py -t $trees -w ${trees.baseName}_output.weights.csv.gz -g guaraniBL -g guaraniNL -g tuca -g iguassuensis --groupsFile $groups
    mv .command.sh .twisst.command.sh
    mv .command.log .twisst.command.log"""
}

process DISTANCE{
    label 'twiss'
    label 'low_resources'
    publishDir "${params.publishDir}/population/phylogeny/02_network/" , mode: 'copy'

    input:
    path (phased_genotypes)

    output:
    path 'genotypes.dist', emit: dmatrix
    tuple path(".network.command.sh"), path(".network.command.log")

    script:
    """
    parseVCF.py -i $phased_genotypes --skipIndels > phased.geno
    distMat.py -g phased.geno -f phased --windType cat -o genotypes.dist -T $task.cpus
    mv .command.sh .network.command.sh
    mv .command.log .network.command.log"""
}

workflow TOPOLOGY_WEIGHTS {
    take:
        GENOTYPES
    main:

    winsize = channel.from(["25","50","100"])
    TREES = ESTIMATE_TREES(GENOTYPES, winsize).trees
    GROUPS = Channel.fromPath(params.groups, checkIfExists: true)
    TWISST(TREES, GROUPS)
    }


workflow PHYLONETWORK{
    take:
        GENOTYPES
    main:
        DISTANCE(GENOTYPES)
}

params.phased = false

include { PHASE_VCF; FIXVCF } from './modules/phasing.nf'  

workflow{
    REFERENCE = Channel.fromPath(params.genome, checkIfExists:true)
    if (params.phased == false){ 
        if (params.phasing == "hifi"){
            VCF = Channel.fromPath(params.vcf, checkIfExists: true)
            BAMS = Channel.fromPath(params.bams, checkIfExists:true).collect()
            VCF = FIXVCF(VCF).vcf
            PHASED_GENO = PHASE_VCF(REFERENCE, VCF, BAMS).phased_vcf
        } else {
            PHASED_GENO = PHASING(VCF).geno
        }
    } else {
        PHASED_GENO = Channel.fromPath(params.phased, checkIfExists:true)
    }
    TOPOLOGY_WEIGHTS(PHASED_GENO)
    PHYLONETWORK(PHASED_GENO)
}