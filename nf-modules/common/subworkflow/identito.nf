/*
 * Identito monitoring
 */

include { identitoPolym } from '../process/identito/identitoPolym'
include { identitoCombine } from '../process/identito/identitoCombine'
include { identitoClustering } from '../process/identito/identitoClustering'

workflow identitoFlow {
    take:
    bam
    fasta
    fai
    polymBed

    main:
    chVersions = Channel.empty()

    identitoPolym(
      bam,
      fasta.collect(),
      fai.collect(),
      polymBed.collect()
    )
    chVersions = chVersions.mix(identitoPolym.out.versions)

    identitoCombine(
      identitoPolym.out.polyms.collect()
    )
    chVersions = chVersions.mix(identitoCombine.out.versions)

    identitoClustering(
      identitoCombine.out.results.collect()
    )

    emit:
    results = identitoClustering.out.results
    versions  = chVersions
}
