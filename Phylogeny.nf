process ESTIMATE_TREES{
    label 'low_resources'
    label 'twiss'
    publishDir "${params.publishDir}/population/phylogeny/01_topologyWeights/" , mode: 'copy'

    input:
    path(phased_genotypes)
    each size
    val(optional)

    output:
    path("phased.geno.phyml_bionj.w*.trees.gz"), emit: trees
    path("phased.geno.phyml_bionj.w*.tsv")
    tuple path(".phyml.command.log"),path(".phyml.command.sh")

    script:
    """
    phyml_sliding_windows.py -T $task.cpus -g $phased_genotypes --prefix phased.geno.phyml_bionj.w${size} -w $size --windType sites $optional --model GTR --optimise n
    mv .command.sh .phyml.command.sh
    mv .command.log .phyml.command.log"""
}

process TWISST {
    label 'medium_resources'
    publishDir "${params.publishDir}/population/phylogeny/01_topologyWeights/" , mode: 'copy'

    label 'twiss'

    input:
    each path(trees)
    each path(groups)
    val (groups_text)
    val (outgroup)

    output:
    path("*.weights.csv.gz")
    tuple path(".*twisst.command.sh"), path(".*twisst.command.log")

    
    script:
    """
    twisst.py -t $trees -w ${trees.baseName}_${outgroup}_output.weights.csv.gz $groups_text --groupsFile $groups --outgroup $outgroup
    mv .command.sh .${trees.baseName}_${outgroup}.twisst.command.sh
    mv .command.log .${trees.baseName}_${outgroup}.twisst.command.log"""
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
    distMat.py -g $phased_genotypes -f phased --windType cat -o genotypes.dist -T $task.cpus
    mv .command.sh .network.command.sh
    mv .command.log .network.command.log"""
}

workflow TOPOLOGY_WEIGHTS {
    take:
        GENOTYPES
        winsize
    main:

    winsize = channel.from(winsize)
    TREES = ESTIMATE_TREES(GENOTYPES, winsize, params.optional).trees
    GROUPS = Channel.fromPath(params.groups, checkIfExists: true)
    TWIST_GROUPS = channel.from(params.twisstgroups)
    OUTGROUP = Channel.from(params.outgroup)
    TWISST(TREES, GROUPS, TWIST_GROUPS, OUTGROUP)
    }


workflow PHYLONETWORK{
    take:
        GENOTYPES
    main:
        DISTANCE(GENOTYPES)
}

params.phased = false
params.twistgroups =  "-g guaraniBL -g guaraniNL -g tuca -g iguassuensis"
params.optional = ""

include { PHASE_VCF; FIXVCF } from './modules/phasing.nf'  

workflow{
    REFERENCE = Channel.fromPath(params.genome, checkIfExists:true)
    if (params.phased == false){ 
        if (params.phasing == "hifi"){
            VCF = Channel.fromPath(params.vcf, checkIfExists: true)
            BAMS = Channel.fromPath(params.bams, checkIfExists:true).collect()
            VCF = FIXVCF(VCF).vcf
            PHASED_GENO = PHASE_VCF(REFERENCE, VCF, BAMS).geno
        } else {
            PHASED_GENO = PHASING(VCF).geno
        }
    } else {
        PHASED_GENO = Channel.fromPath(params.phased, checkIfExists:true)
    }
    TOPOLOGY_WEIGHTS(PHASED_GENO, ["25","50","100"])
    PHYLONETWORK(PHASED_GENO)
}