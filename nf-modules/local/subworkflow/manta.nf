/*
 * Manta structural variant calling Flow
 */

include { manta } from '../../local/process/manta'

workflow mantaFlow {

    take:
    bamTN
    bam
    bed
    fasta
    fastaFai

    main:
    chVersions = Channel.empty()

    // Mutect2 tumor/bam inputs
    //[meta][tumor_bam, normal_bam],[tumor_bai, normal_bai]
    bamTN
      .map{ it -> [it[0], [it[1], it[3]], [it[2], it[4]]] }
      .set{ chBamManta }

    chBamMantaCombined = chBamManta.mix(bam)
    chBamMantaCombined.view()

      manta(
        chBamMantaCombined,
        bed,
        fasta,
        fastaFai
        )

    emit:
    versions = chVersions
}
