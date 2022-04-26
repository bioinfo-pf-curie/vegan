/*
 * QC Worflow for Bam files
 */

 include { getFragmentSize } from '../../common/process/getFragmentSize'
 include { getSeqDepth } from '../../common/process/getSeqDepth'
 include { prepareExonInfo } from '../../local/process/prepareExonInfo'
 include { genesCoverage } from '../../common/process/genesCoverage'
 include { getWGSmetrics } from '../../common/process/getWGSmetrics'

workflow bamQcFlow2 {
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

    getSeqDepth.out.versions.view()

    prepareExonInfo(
      bed.collect(),
      gtf.collect()
      )

    genesCoverage(
      bamFiltered,
      prepareExonInfo.out.exonBed
      )

    getWGSmetrics(
      bamFiltered,
      bed,
      fasta,
      dict
      )

    chVersions = chVersions.mix(getSeqDepth.out.versions)
    chVersions.view()

  emit:
  fragSize = getFragmentSize.out.metrics
  seqDepth = getSeqDepth.out.metrics
  bedDepth = getSeqDepth.out.mosdepthBed
  versions = getFragmentSize.out.versions
}
