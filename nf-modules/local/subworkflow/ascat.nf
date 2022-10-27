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

    countsCh
      .combine(countsCh)
      .filter { it[0].pair_id == it[2].pair_id && it[0].status == "tumor" && it[2].status == "normal"}
      .map{ it ->
        meta = [tumor_id:it[0].id, normal_id:it[2].id, status: "pair", id:it[2].pair_id, sex:it[0].sex]
      return [meta, it[1], it[3]]}
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
