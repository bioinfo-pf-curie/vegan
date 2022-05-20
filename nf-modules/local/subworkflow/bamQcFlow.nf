/*
 * QC Worflow for Bam files
 */

include { getFragmentSize } from '../../common/process/getFragmentSize'
include { getSeqDepth } from '../../common/process/getSeqDepth'
include { prepareExonInfo } from '../../local/process/prepareExonInfo'
include { genesCoverage } from '../../common/process/genesCoverage'
include { getWGSmetrics } from '../../common/process/getWGSmetrics'

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
      getFragmentSize(
        bamFiltered)

      chVersions = chVersions.mix(getFragmentSize.out.versions)
    }

    getSeqDepth(
      bamFiltered,
      bed.collect()
    )

    prepareExonInfo(
      bed.collect(),
      gtf.collect()
    )

    genesCoverage(
      bamFiltered,

      prepareExonInfo.out.exonBed.collect()
      )

    getWGSmetrics(
      bamFiltered,
      bed.collect(),
      fasta.collect(),
      dict.collect( )
      )

    chVersions = chVersions.mix(getSeqDepth.out.versions)

  emit:
    fragSize = getFragmentSize.out.metrics
    seqDepth = getSeqDepth.out.metrics
    bedDepth = getSeqDepth.out.mosdepthBed
    versions = getFragmentSize.out.versions
}
