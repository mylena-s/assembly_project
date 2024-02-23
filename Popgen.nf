
log.info """\
    P I P E L I N E     F O R    P O P G E N  
    ==========================================
   
    performs long read mapping with minimap
    sorts, marks duplicates, and calls variants    
    Current project dir. is: $projectDir"""
    .stripIndent()

include { MINIMAP; BAMQC; FREEBAYES_JOINT; MARKDUP } from './modules/minimap.nf'  
params.dtype = "hifi"
params.mapped = false
params.nodup = false
workflow {
    READS = Channel.fromPath(params.reads, checkIfExists:true) 
    ASSEMBLY = Channel.fromPath(params.genome, checkIfExists: true)    
    INPUT = ASSEMBLY.combine(READS)  
    COLLECTED_MAPPINGS = channel.empty()
    if (params.mapped == false) {
        MAPPINGS = MINIMAP(INPUT).mappings
    } else {
        MAPPINGS = Channel.fromPath(params.mapped, checkIfExists: true) 
    }
    if (params.nodup == false) {
        MAPPINGS_NODUP = MARKDUP(MAPPINGS).mappings
        COLLECTED_MAPPINGS = MAPPINGS_NODUP.collect()
    } else {
        MAPPINGS_NODUP = Channel.fromPath(params.nodup, checkIfExists:true) }
        COLLECTED_MAPPINGS = MAPPINGS_NODUP.collect() 
    FREEBAYES_JOINT(ASSEMBLY, COLLECTED_MAPPINGS)
}


