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
  targetBed
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
    chFragSize = collectInsertSizeMetrics.out.results
  }else{
    chFragSize = Channel.empty()
  }

  mosdepth(
    bamFiltered,
    targetBed
  )
  chVersions = chVersions.mix(mosdepth.out.versions)
  if (params.targetBed){
    chMosdepthLog = mosdepth.out.regionsTxt
  }else{
    chMosdepthLog = mosdepth.out.globalTxt
  }

  prepareExonInfo(
    gtf,
    targetBed
  )

  genesCoverage(
    bamFiltered,
    prepareExonInfo.out.exonBed.collect()
  )

  collectWgsMetrics(
    bamFiltered,
    targetBed,
    fasta,
    dict
  )
  chVersions = chVersions.mix(collectWgsMetrics.out.versions)

  emit:
    fragSize = chFragSize
    depth = chMosdepthLog
    geneCovMqc = genesCoverage.out.geneCovMqc
    wgsMetrics = collectWgsMetrics.out.metrics
    versions = chVersions
}
