#!/usr/bin/env nextflow

log.info """\
    P I P E L I N E    F O R   R E P E A T 
              A N N O T A T I O N 
    =======================================
    Uses NextDenovo as ONT assembler or Hifiasm
                 for Pacbio hifi. 
    
    Also performs QC with different tools
    
    Current project dir. is: $projectDir"""
    .stripIndent()

process EARLGREY {
    label 'earlGrey'
    label 'resource_intensive'
    publishDir "${params.publishDir}/annotation/repeat/", mode: 'copy'
    
    input:
    tuple path(genome), path(sat)

    output:
    path '02_earlGrey_*'
    path "02_earlGrey_*/${genome.baseName}.masked", emit: masked
    path "02_earlGrey_*/${genome.baseName}.library", emit: library

    script:
    """
    cp $genome ${genome.baseName}_copy.fasta
    sed -E 's/>(.*)#(.*)\\s.*/>\\2#Satellite\\/\\2/g' $sat | sed 's/_reads//g' > satellites.fasta
    earlGrey -l \$PWD'/satellites.fasta' -g ${genome.baseName}_copy.fasta -s ${genome.baseName} -o ./02_earlGrey -t $task.cpus -d yes
    mv .command.log .command.sh 02_earlGrey/${genome.baseName}
    """
}  

process MASK {
    label 'earlGrey'
    label 'resource_intensive'
    publishDir "${params.publishDir}/annotation/repeat/", mode: 'copy'

    input:
    tuple path(genome), path(lib)

    output:
    path "04_masking"
    path "04_masking/*masked", emit: masked_genome
    script:
    """
    RepeatMasker $genome -lib $lib -xsmall -norna -s -pa $task.cpus -a -dir 04_masking
    calcDivergenceFromAlign.pl -s ${genome.simpleName}.divsum 04_masking/*.align
    mv *divsum 04_masking
    mv .command.sh .command.log 04_masking"""
}
    
process BLAST_TES{
    label 'repeats'
    label 'medium_resources'

    input:
    path(library)
    output:
    path("${library.simpleName}_pfam.out")

    script:
    """
    getorf -sequence $library -outseq out.orf -minsize 300
    pfam_scan.pl -fasta out.orf -dir $projectDir/lib/Pfam/ -outfile ${library.simpleName}_pfam.out
    """

}

process SAT_ANNOT {
    label 'minimap'
    label 'resource_intensive'
    publishDir "${params.publishDir}/annotation/repeat/", mode: 'copy'

    input:
    path genome
    output:
    path "01_srf_${genome.baseName}"
    path "01_srf_${genome.baseName}/*.srf.fa", emit: lib
    path "01_srf_${genome.baseName}/*bed", emit: satbed

    script:
    """
    mkdir -p tmp_dir
    kmc -fm -k151 -t16 -ci20 -cs100000 $genome count.kmc tmp_dir
    kmc_dump count.kmc count.txt
    srf -p prefix count.txt > ${genome.baseName}.srf.fa
    minimap2 -c -N1000000 -f1000 -r100,100 <(srfutils.js enlong ${genome.baseName}.srf.fa) $genome > srf-aln.paf
    srfutils.js paf2bed srf-aln.paf > ${genome.baseName}.srf-aln.bed
    srfutils.js bed2abun -g $params.gs ${genome.baseName}.srf-aln.bed > ${genome.baseName}.srf-aln.len
    mkdir 01_srf_${genome.baseName}
    cp .command.sh 01_srf_${genome.baseName}/.${genome.baseName}.command.sh
    cp .command.log 01_srf_${genome.baseName}/.${genome.baseName}.command.logs
    mv *.bed *.len 01_srf_${genome.baseName}/
    mv ${genome.baseName}.srf.fa 01_srf_${genome.baseName}/
    """
}

process CENTROMERE{
    label 'medium_resources'
    publishDir "${params.publishDir}/annotation/repeat/", mode: 'copy', pattern: '03_nessie*'
    errorStrategy 'ignore'
    label 'minimap'
    input:
    path genome
    each mode

    output:
    path '*_complexity.tsv', optional: true, emit: complexity
    path '*_entropy.tsv', optional: true, emit: entropy

    script:
    """
    mkdir 03_nessie_${genome.baseName}
    nessie -I $genome -O $genome$mode-l 10000 -s 1000
    for file in *tsv; do to_wig.py -i \$file -o \$file'.wig'; done
    for file in *wig; do tail -n +2 \$file | grep -v 'track type' - > \$file'.bw'; done
    samtools faidx $genome
    awk '{print \$1, \$2}' $genome'.fai' > $genome'.chrsizes'
    for file in *.bw; do wigToBigWig \$file *.chrsizes \$file'.bigwig'; done
    cp *tsv 03_nessie_${genome.baseName}
    cp *bigwig 03_nessie_${genome.baseName}
    cp .command.sh 03_nessie_${genome.baseName}/.${genome.baseName}.command.sh
    cp .command.log 03_nessie_${genome.baseName}/.${genome.baseName}.command.log
    """
}

process BRAKER2 { 
    label 'braker'
    label 'resource_intensive'
    publishDir "${params.publishDir}/annotation/gene/", mode: 'copy'
    
    input:
    path genome
    path proteins

    output:
    
    script:
    """
    gunzip -c $proteins > proteins.fa
    braker.pl --species=Crenicichla --genome=$genome --prot_seq=proteins.fa --workingdir=braker2_${genome.simpleName} --threads=$task.cpus --busco_lineage=actinopterygii_odb10
    """
}
params.gs=800000
params.satlib = false
params.centromere = false
params.replib = false
params.masked = false
params.proteins = "lib/Actinopterygii.fasta.gz" 

workflow GENE_ANNOTATION{
    take:
        GENOME
    
    main:
        REFPROT = Channel.fromPath(params.proteins,  checkIfExists: true)
        INPUT = ASSEMBLIES.combine(REFPROT)
        ANNOT = BRAKER2(INPUT)
}

workflow REPEAT_ANNOTATION{
    take:
        GENOME
    main:
        if (params.satlib == false){
            SATELLITE = SAT_ANNOT(GENOME).lib
        }else{
            SATELLITE = Channel.fromPath(params.satlib, checkIfExists:true)
        }
        if (params.replib==false){
            INPUT = GENOME.combine(SATELLITE)
            REPEATS = EARLGREY(INPUT)
            MASKED = REPEATS.masked
            LIBRARY = REPEATS.library
            BLAST_TES(LIBRARY)

        } else {
            LIBRARY = Channel.fromPath(params.replib)
            INPUT = GENOME.combine(LIBRARY)
            MASKED = MASK(INPUT)
        }
        if (params.centromere) {
            modes = ['_complexity.tsv -L ','_entropy.tsv -E ']
            CENTROMERE(GENOME, modes)}
    emit:
        MASKED
}

workflow {
    ASSEMBLIES = Channel.fromPath(params.genome, checkIfExists: true)
    
    if (params.masked == false) {
        MASKED = REPEAT_ANNOTATION(ASSEMBLIES)
    } else {
        MASKED = ASSEMBLIES
    } 
    
    GENE_ANNOTATION(MASKED)
}
