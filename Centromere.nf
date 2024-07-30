include { CENTROMERE } from './modules/centromere.nf'

workflow {
    ASSEMBLIES = Channel.fromPath(params.genome, checkIfExists: true)
    modes = ['_complexity.tsv -L ','_entropy.tsv -E ']
    CENTROMERE(ASSEMBLIES, modes)
}
