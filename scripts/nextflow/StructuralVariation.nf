include { MINIMAP2; BAMQC} from './modules/minimap.nf'  

process DEBREAK {
    publishDir "${params.publishDir}/structural/debreak", mode: 'copy'
    label 'structural'
    label 'resource_intensive'
    input:
    path(genome)
    each path(bam)

    output:
    path("*debreak_out")
    path("*debreak_out/debreak.vcf"), emit: svs

    script:
    """
    samtools index $bam
    debreak --bam $bam -d $params.covmax -m $params.covmin -o ${bam.simpleName}_debreak_out/ --rescue_large_ins --rescue_dup --poa --ref $genome -t $task.cpus
    cp .command.sh .command.log ${bam.simpleName}_debreak_out/
    """
}

process VACMAP {
    publishDir "${params.publishDir}/structural/alignment", mode: 'copy'
    label 'vacmap'
    label 'resource_intensive'
    
    input:
    tuple path(genome), path(reads)

    output:
    path("*sorted.bam"), emit: mappings
    tuple path(".vacmap.command.sh"), path(".vacmap.command.log")

    script:
    """
    vacmap -ref $genome -read $reads -mode $params.mode -t $task.cpus -workdir . --H --fakecigar > alignment.sam
    samtools sort -@$task.cpus alignment.sam > ${reads.simpleName}.sorted.bam
    mv .command.sh .vacmap.command.sh
    mv .command.log .vacmap.command.log
    """

}

process SNIFFLES{
    publishDir "${params.publishDir}/structural/sniffles", mode: 'copy'
    label 'sniffles'
    label 'resource_intensive'
    
    input:
    path(genome)
    each path(bam)

    output:
    path("${bam.simpleName}_sniffles_out")
    path("${bam.simpleName}_sniffles_out/sniffles.vcf"), emit: svs
    val("${bam.simpleName}"), emit: name

    script:
    """
    samtools index $bam
    mkdir ${bam.simpleName}_sniffles_out
    sniffles --input $bam --vcf sniffles.vcf --threads $task.cpus --reference $genome
    mv sniffles.vcf ${bam.simpleName}_sniffles_out
    cp .command.sh .command.log ${bam.simpleName}_sniffles_out/
    """
}

process SNIFFLES_MULTI{
    publishDir "${params.publishDir}/structural/sniffles", mode: 'copy'
    label 'sniffles'
    label 'resource_intensive'
    
    input:
    path(genome)
    each path(bam)

    output:
    path("*snf"), emit: snif
    tuple path(".${bam.simpleName}.command.sh"), path(".${bam.simpleName}.command.log"), emit: logs

    script:
    """
    samtools index $bam
    mkdir ${bam.simpleName}_sniffles_out
    sniffles --input $bam --snf ${bam.simpleName}.snf --threads $task.cpus --reference $genome
    cp .command.sh .${bam.simpleName}.command.sh
    cp .command.log .${bam.simpleName}.command.log
    """
}


process SNIFFLES_JOINTCALL{
    publishDir "${params.publishDir}/structural/sniffles", mode: 'copy'
    label 'sniffles'
    label 'resource_intensive'
    
    input:
    path(sniff)

    output:
    path("*vcf"), emit: vcf
    tuple path(".sniffles.command.sh"), path(".sniffles.command.log"), emit: logs

    script:
    def files = sniff.join(' ')
    """
    sniffles --input $files --threads $task.cpus --vcf structural_jointcall.vcf
    cp .command.sh .sniffles.command.sh
    cp .command.log .sniffles.command.log
    """
}

process WINNOWMAP {
    publishDir "${params.publishDir}/structural/alignment", mode: 'copy'
    label 'minimap'
    label 'resource_intensive'
    
    input:
    tuple path(genome), path(reads)
    val(type)

    output:
    path("*sorted.bam"), emit: mappings
    tuple path(".winnowmap.command.sh"), path(".winnowmap.command.log")

    script:    
    """meryl count threads=$task.cpus k=15 output merylDB $genome
meryl print greater-than distinct
=0.9998 merylDB > repetitive_k15.txt
winnowmap -t $task.cpus -W repetitive_k15.txt -ax map-$type $genome $reads > output.sam
samtools view -q 15 -S -b output.sam > unsorted.bam
samtools sort -@${task.cpus} -O BAM -o ${genome.baseName}_${reads.simpleName}_winnowmap_sorted.bam unsorted.bam
rm unsorted.bam output.sam 
cp .command.sh .winnowmap.command.sh
cp .command.log .winnowmap.command.log
"""}


process MERGE{
    publishDir "${params.publishDir}/structural/", mode: 'copy'
    label 'jasmine'
    label 'low_resources'

    input:
    tuple path(SV1), path(SV2)
    path(genome)
    val(sample)

    output:
    path ("*_merged.vcf"), emit: merged
    tuple path(".*command.sh"), path(".*command.log")
    script:
    """jasmine --threads $task.cpus --normalize_type --comma_filelist file_list=$SV1,$SV2 out_file=${sample}_merged.vcf --ignore_strand --dup_to_ins genome_file=$genome
    cp .command.sh .jasmine.command.sh
    cp .command.log .jasmine.command.log"""
}

process SUMMARY {
    publishDir "${params.publishDir}/structural/alignment", mode: 'copy'
    input:
    path(merged_vcf)
    output:
    path("merged_summary_table.csv"), emit: summary
    script:
    """#!/usr/bin/env python
import csv
from collections import defaultdict

# Initialize counts
counts = defaultdict(lambda: defaultdict(int))

# Input and output file names
input_file = '${merged_vcf}'
output_file = "merged_summary_table.csv"

# Process the file
with open(input_file, "r") as file:
    for line in file:
        # Skip empty or malformed lines
        if not line.strip():
            continue
        
        # Parse the line
        columns = line.strip().split("\\t")
        if len(columns) < 8:  # Ensure the line has enough columns
            continue
        
        # Extract relevant fields
        info_field = columns[7]  # INFO field is in the 8th column (index 7)
        info_data = {key: value for key, value in (item.split("=") for item in info_field.split(";") if "=" in item)}
        
        # Get SVTYPE and SUPP_VEC
        svtype = info_data.get("SVTYPE", "NA")
        supp_vec = info_data.get("SUPP_VEC", "NA")
        # Count occurrences if SUPP_VEC is relevant
        if supp_vec in {"10","01", "11"}:
            counts[svtype][supp_vec] += 1

# Write results to a CSV file
with open(output_file, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    # Write header
    writer.writerow(["SVTYPE", "SUPP_VEC=10","SUPP_VEC=01", "SUPP_VEC=11"])
    for svtype in sorted(counts.keys()):
        writer.writerow([svtype, counts[svtype]["10"], counts[svtype]["01"], counts[svtype]["11"]])
"""
}

// def findCoverage {
//     input:
//     path txt

//     output:
//     val(coverage), emit coverage

//     script:
//     """value=\$(grep -oP '(?<=mean coverageData = )\\d+\\.\\d+' $txt)
// echo \$value
//     """
// }

params.mapped = false
params.mode = "L"
params.type = "ont"
params.sniffles = false
params.debreak = false
params.multi = false

workflow {
    ASSEMBLY = Channel.fromPath(params.genome)    
    if (params.mapped == false) {
        READS = Channel.fromPath(params.reads, checkIfExists:true)
        INPUT = ASSEMBLY.combine(READS)
    //    MAPPINGS = VACMAP(INPUT).mappings // vacmap is giving weird results
        // MAPPINGS = WINNOWMAP(INPUT, params.type).mappings
        MAPPINGS = MINIMAP2(INPUT, params.type).mappings
    } else {
        MAPPINGS = Channel.fromPath(params.mapped)
    }
    // txt = BAMQC(MAPPINGS).txt
    // COVERAGE = findCoverage(txt).coverage
    if (params.multi == false) {
        if (params.debreak != false){
            SVS1 = DEBREAK(ASSEMBLY, MAPPINGS).svs}
        if (params.sniffles != false){
            SVS2 = SNIFFLES(ASSEMBLY, MAPPINGS)
            if (params.debreak != false){ 
                SVSS = SVS1.combine(SVS2.svs)
                MERGEDSVS = MERGE(SVSS, ASSEMBLY, SVS2.name) 
                TXT = SUMMARY(MERGEDSVS.merged)
            }}
    }else{
        SNIFS = SNIFFLES_MULTI(ASSEMBLY, MAPPINGS).snif.collect()
        SVSS = SNIFFLES_JOINTCALL(SNIFS)
    }
}
