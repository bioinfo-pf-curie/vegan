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

    // Split the paired BAM
    // tumor_id, normal_id, pair_id, id, status:pair, sex
    chTumorBam = bam.map{ it ->
        def meta = [id:it[0].tumor_id, pair_id:it[0].pair_id, status: "tumor", sex:it[0].sex]
        return [meta, it[1], it[2] ]
      }

    chNormalBam = bam.map{ it ->
        def meta = [id:it[0].normal_id, pair_id:it[0].pair_id, status: "normal", sex:it[0].sex]
        return [meta, it[3], it[4] ]
      }

    alleleCounter(
      chNormalBam.mix(chTumorBam),
      acLoci,
      dict,
      fasta,
      fastaFai,
    )

    countsCh = alleleCounter.out.counts
      .combine(alleleCounter.out.counts)
      .filter { it[0].pair_id == it[2].pair_id && it[0].status == "tumor" && it[2].status == "normal"}
      .map{ it ->
        meta = [tumor_id:it[0].id, normal_id:it[2].id, status: "pair", pair_id:it[2].pair_id, id:it[0].id+"_"+it[2].id, sex:it[0].sex]
        return [meta, it[1], it[3]]
      }

    convertAlleleCounts(
      countsCh
    )

    ascat(
      convertAlleleCounts.out.counts,
      acLociGC
    )

    emit:
    results = ascat.out.results
    versions = chVersions
}
