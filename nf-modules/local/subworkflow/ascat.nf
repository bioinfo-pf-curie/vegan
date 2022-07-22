/*
 * Ascat flow
 */

include { alleleCounter } from '../../common/process/ascat/alleleCounter'
include { convertAlleleCounts } from '../../common/process/ascat/convertAlleleCounts'
include { ascat } from '../../common/process/ascat/ascat'

workflow ascatFlow {

    take:
    bam
    acLoci
    acLociGC
    dict
    fasta
    fastaFai

    main:
    chVersions = Channel.empty()

    alleleCounter(
      bam,
      acLoci,
      dict,
      fasta,
      fastaFai,
      )

    countsCh = alleleCounter.out.alleleCounterOutCh

    countsCh.combine(countsCh)
    .filter {it[0].tumor_id == it[1].baseName && it[0].normal_id == it[3].baseName}
    .map { it ->  return [it[0], it[1], it[3] ]}
    .set{ countsCombinedCh }

    convertAlleleCounts(
      countsCombinedCh
      )

    ascat(
      convertAlleleCounts.out.counts,
      acLociGC
      )

    emit:

    versions = chVersions
}
