process MINIMAP2{
    label 'resource_intensive'
    label 'syri'
    maxForks 2
    publishDir "${params.publishDir}/alignment", mode: 'copy' 

    input:
    tuple path(genome), val(number)
    
    output:
    tuple path ('*.sam'), path("${genome[0]}", includeInputs: true), path("${genome[1]}",includeInputs: true), val(number),  emit: mappings
    path '.*command.*'

    script:
    def name = "${genome[0].simpleName}_${genome[1].simpleName}"

    """minimap2 -ax asm20 -t $task.cpus --eqx ${genome[0]} ${genome[1]} > ${name}.sam
    mv .command.sh .${name}.command.sh
    mv .command.log .${name}.command.log
    """
}

process SYRI{
    label 'medium_resources'
    label 'syri'
    maxForks 10
    publishDir "${params.publishDir}/synteny", mode: 'copy'   

    input:
    tuple path(sam), path(reference), path(query), val(number)

    output:
    tuple path(".*command.sh"), path(".*command.log")
    tuple path("*syri.out"), val(number), emit: variants
    path ("*summary")

    script:
    """syri -c $sam -r $reference -q $query -k -F S -f --nosnp --prefix ${sam.simpleName}"""


}
process SYNTENY{
    label 'medium_resources'
    label 'syri'
    maxForks 10
    publishDir "${params.publishDir}/synteny", mode: 'copy'   

    input:
    tuple path(syri1), path(syri2)
    path("genomes.txt")

    output:
    tuple path(".*command.sh"), path(".*command.log")
    path("*svg")

    script:
    """plotsr --itx --sr $syri1 --sr $syri2 --genomes genomes.txt --chr LG1 --chr LG2 --chr LG3 --chr LG4 --chr LG5 --chr LG6 -H 3 -W 10 -s 5000 -S 0.5 -f 12 -R --cfg $projectDir'/config.cfg' --tracks $projectDir'/tracks.txt' -o synteny_part1.svg -b svg
    plotsr --itx --sr $syri1 --sr $syri2 --genomes genomes.txt --chr LG8 --chr LG9 --chr LG10 --chr LG11 --chr LG12 --chr LG13 -H 3 -W 10 -s 5000 -S 0.5 -f 12 -R --cfg $projectDir'/config.cfg' --tracks $projectDir'/tracks.txt' -o synteny_part2.svg -b svg
    plotsr --itx --sr $syri1 --sr $syri2 --genomes genomes.txt --chr LG14 --chr LG16 --chr LG18 --chr LG19 --chr LG20 --chr LG22 -H 3 -W 10 -s 5000 -S 0.5 -f 12 -R --cfg $projectDir'/config.cfg' --tracks $projectDir'/tracks.txt' -o synteny_part3.svg -b svg
    plotsr --itx  --sr $syri1 --sr $syri2 --genomes genomes.txt --chr LG7 --chr LG21 --chr LG15 --chr LG17 --chr LG23 --chr LG24 -H 3 -W 10 -s 5000 -S 0.5 -f 12 -R --cfg $projectDir'/config.cfg' --tracks $projectDir'/tracks.txt' --markers $projectDir'/markers.bed' -o synteny_part4.svg -b svg
    """
}
workflow {
    reference = Channel.fromPath(params.reference) 
    query1 = Channel.fromPath(params.query)
    query2 = Channel.fromPath(params.query2)
    genomes = Channel.fromPath(params.genomes_file)
    
    // Create genome pairs with numbers
    genome_files = reference.concat(query1, query2)
                .collate(2, 1, remainder= false)
                .merge(channel.from(1, 2)) { a, b -> tuple(a, b) }
    
    // Run alignment
    alignments = MINIMAP2(genome_files)
    
    // Run SYRI and collect outputs
    syri_results = SYRI(alignments.mappings)
    
    // Properly pair the SYRI outputs for SYNTENY
    // Assuming you want to compare reference vs query1 and reference vs query2
    syri_pairs = syri_results.variants
        .toSortedList { a, b -> a[1] <=> b[1] }
        .map { list -> 
            // This assumes the first output is reference vs query1 (number=1)
            // and second is reference vs query2 (number=2)
            tuple(list[0][0], list[1][0]) 
        }
    
    // Run SYNTENY with the paired SYRI outputs
    SYNTENY(syri_pairs, genomes)
}
