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
    runName

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
      identitoPolym.out.tsv.collect(),
      runName
    )
    chVersions = chVersions.mix(identitoCombine.out.versions)

    identitoClustering(
      identitoCombine.out.tsv.collect(),
      runName
    )

    emit:
    tsv = identitoCombine.out.tsv
    results = identitoClustering.out.results
    versions  = chVersions
}
