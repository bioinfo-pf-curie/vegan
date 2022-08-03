/*
 * Manta structural variant calling Flow
 */

include { manta } from '../../local/process/manta'

workflow mantaFlow {

    take:
    bamTN
    bed
    fasta
    fastaFai

    main:
    chVersions = Channel.empty()

      manta(
        bamTN,
        bed,
        fasta,
        fastaFai
        )

    emit:
    versions = chVersions
}
