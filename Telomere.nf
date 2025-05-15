process FIND_TELOMERES{
    label 'medium_resources'
    label 'tel'
    maxForks 2
    publishDir "${params.publishDir}/telomere", mode: 'copy'   

    input:
    path(genome)

    output:
    path("*_telomeric_repeat_windows.tsv"), emit: tel
    tuple path(".*command.sh"), path(".*command.log")

    script:
    def string = "TTAGGG" 
"""seqkit seq -u $genome > upper.genome.fasta
tidk search -s $string -w 1000 -o ${genome.simpleName}_${string} --dir . upper.genome.fasta
mv .command.sh .telomere.command.sh
mv .command.log .telomere.command.log
"""
}

process tsv2Bed{
    label 'syri'
    publishDir "${params.publishDir}/telomere", mode: 'copy'   

    input:
    path(tel)

    output:
    path("*bedgraph")
    tuple path(".*command.sh"), path(".*command.log")

    script:
    """#!/usr/bin/env python3
import pandas as pd
df = pd.read_csv("$tel", sep="\\t")
df["start"] = df.window-10000
df["count"] = df.forward_repeat_number + df.reverse_repeat_number
df.loc[:,["id","start","window","count"]].to_csv("${tel.simpleName}.bedgraph", sep="\\t", header=None, index=None)"""
}


workflow {
    genome = Channel.fromPath(params.reference)
    telomeres = FIND_TELOMERES(genome).tel
    tsv2Bed(telomeres)
    }
