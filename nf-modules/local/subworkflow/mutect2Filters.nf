/*
 * Mutect2 Filters
 */

include { getPileupSummaries } from '../../common/process/gatk/getPileupSummaries'
include { gatherPileupSummaries } from '../../common/process/gatk/gatherPileupSummaries'
include { calculateContamination } from '../../common/process/gatk/calculateContamination'
include { filterMutect2Calls } from '../../common/process/gatk/filterMutect2Calls'

workflow mutect2FiltersFlow {

  take:
  bam
  intervals
  germlineResource
  germlineResourceIndex
  bed
  dict
  mutect2CallsToFilter
  stats
  fasta
  fastaFai

  main:
  chVersions = Channel.empty()

  getPileupSummaries(
    bam,
    intervals,
    germlineResource,
    germlineResourceIndex,
    bed
    )

  gatherPileupSummaries(
    getPileupSummaries.out.pileupSummaries,
    dict
    )

  pairBamCalculateContaminationCh = bam.join(gatherPileupSummaries.out.mergedPileupFileCh)

  calculateContamination(
    pairBamCalculateContaminationCh
    )

  mutect2CallsToFilter = mutect2CallsToFilter.map{
      meta, variantCaller, vcf, index ->
      [meta, vcf, index]
  }.join(stats)

  if (!params.skipMutectContamination){
    mutect2CallsToFilter = mutect2CallsToFilter.join(calculateContamination.out.contaminationTable)
  }else{
    mutect2CallsToFilter = mutect2CallsToFilter.combine(Channel.from('NO_FILE'))
  }

  //mutect2CallsToFilter.view()

  filterMutect2Calls(
    mutect2CallsToFilter,
    dict,
    fasta,
    fastaFai,
    germlineResource,
    germlineResourceIndex,
    intervals
    )



  emit:
  versions = chVersions
}
