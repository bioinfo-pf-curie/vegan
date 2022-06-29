/*
 * QC Worflow for Bam files
 */

include { collectInsertSizeMetrics } from '../../common/process/picard/collectInsertSizeMetrics'
include { mosdepth } from '../../common/process/mosdepth/mosdepth'
include { prepareExonInfo } from '../../local/process/prepareExonInfo'
include { genesCoverage } from '../../local/process/genesCoverage'
include { collectWgsMetrics } from '../../common/process/gatk/collectWgsMetrics'

workflow bamQcFlow {

  take:
    bamFiltered
    bed
    gtf
    fasta
    dict

  main:
    chVersions = Channel.empty()

    if (!params.singleEnd){
      collectInsertSizeMetrics(
        bamFiltered
      )
      chVersions = chVersions.mix(collectInsertSizeMetrics.out.versions)
    }

    mosdepth(
      bamFiltered,
      bed.collect()
    )
    chVersions = chVersions.mix(mosdepth.out.versions)

    prepareExonInfo(
      bed.collect(),
      gtf.collect()
    )

    genesCoverage(
      bamFiltered,
      prepareExonInfo.out.exonBed.collect()
    )

    collectWgsMetrics(
      bamFiltered,
      bed.collect(),
      fasta.collect(),
      dict.collect( )
    )
    chVersions = chVersions.mix(collectWgsMetrics.out.versions)

  emit:
    fragSize = collectInsertSizeMetrics.out.results
    seqDepth = mosdepth.out.metrics
    bedDepth = mosdepth.out.mosdepthBed
    versions = chVersions
}
