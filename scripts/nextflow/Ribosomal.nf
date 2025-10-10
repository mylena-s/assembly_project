process ribotin{
    publishDir "${params.publishDir}/rRNA/" , mode: 'copy'
    label 'ribotin'
    input:
    tuple path(tangle), path(hifiasm_folder), val(hifiasm_prefix)
    tuple path(hifi), path(ont)
    val (size)
    path(ref)

    output:
    path('ribotin_output*'), emit:output
    script:
    """ribotin-hifiasm -a $hifiasm_folder/$hifiasm_prefix -o ribotin_output -i $hifi --nano $ont -c $tangle --approx-morphsize $size --orient-by-reference $ref  
    cp .command* ribotin_output*"""
}

process fastqfasta{
    input:
    tuple path(hifi), path(ont)
    output:
    tuple path("${hifi.simpleName}.fasta"), path("${ont.simpleName}.fasta"), emit:reads

    script:
    """seqtk seq -a $hifi > ${hifi.simpleName}'.fasta'
    seqtk seq -a $ont > ${ont.simpleName}'.fasta'
    """
}

process makepangenome {
    publishDir "${params.publishDir}/rRNA/" , mode: 'copy'

    label 'pggb'
    label 'medium_resources'

    input:
    path(ribotin_output)

    output:
    path("rna_pangenome")

    script:
    """
    cat $ribotin_output/consensus.fa $ribotin_output/morphs.fa > rnacluster.fa
    NSEQ=\$(grep -c ">" rnacluster.fa)
    samtools faidx rnacluster.fa
    pggb -i rnacluster.fa -o rna_pangenome -n \$NSEQ -t $task.cpus -p 90 -s 2000
    cp .command* rna_pangenome"""
}

params.fasta = false
params.morphsize = 13000
workflow{
    HIFIASM = Channel.fromPath(params.hifiasm_folder)
    HIFIASM_P = Channel.value(params.hifiasm_prefix)
    TANGLE = Channel.fromPath(params.tangle)
    hifi = Channel.fromPath(params.hifi)
    ont = Channel.fromPath(params.ont)
    if (params.fasta == false){
        reads = fastqfasta(reads).reads}
    else {
        reads = hifi.combine(ont)
    }
    INPUT = TANGLE.combine(HIFIASM).combine(HIFIASM_P)
    REFERENCE = Channel.fromPath(params.ref)
    RIBOSOME =  ribotin(INPUT, reads, params.morphsize, REFERENCE).output
    makepangenome(RIBOSOME)
}