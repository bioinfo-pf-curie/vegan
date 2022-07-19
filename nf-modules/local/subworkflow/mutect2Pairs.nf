/*
 * Mutect2 Tumor/Normal
 */

include { mutect2 } from '../../common/process/gatk/mutect2'
include { mergeMutect2Stats } from '../../common/process/gatk/mergeMutect2Stats'
include { concatVCF } from '../../local/process/concatVCF'
include { getPileupSummaries } from '../../common/process/gatk/getPileupSummaries'
include { gatherPileupSummaries } from '../../common/process/gatk/gatherPileupSummaries'
include { calculateContamination } from '../../common/process/gatk/calculateContamination'
include { filterMutect2Calls } from '../../common/process/gatk/filterMutect2Calls'

workflow mutect2PairsFlow {

  take:
  bam
  bed
  fasta
  fai
  dict
  germlineResource
  germlineResourceIndex
  panelsOfNormals
  panelsOfNormalsIndex
  intervals

  main:
  chVersions = Channel.empty()

  // Mutect2 tumor/bam inputs
  //[meta][tumor_bam, normal_bam],[tumor_bai, normal_bai]
  bam
    .map{ it -> [it[0], [it[1], it[3]], [it[2], it[4]]] }
    .set{ chBamMutect2 }

  /*
   * MUTECT2 CALLS
   */

  mutect2(
    chBamMutect2,
    bed,
    fasta,
    fai,
    dict,
    germlineResource,
    germlineResourceIndex,
    panelsOfNormals,
    panelsOfNormalsIndex
  )
  chVersions = chVersions.mix(mutect2.out.versions)

  mergeMutect2Stats(
    mutect2.out.stats
  )
  chVersions = chVersions.mix(mergeMutect2Stats.out.versions)

  concatVCF(
    mutect2.out.vcf,
    bed,
    fasta,
    fai
  )
  chVersions = chVersions.mix(concatVCF.out.versions)

  /*
   * PILEUP SUMMARY
   */

  getPileupSummaries(
    bam,
    intervals,
    germlineResource,
    germlineResourceIndex,
    bed
  )
  chVersions = chVersions.mix(getPileupSummaries.out.versions)

  gatherPileupSummaries(
    getPileupSummaries.out.pileupSummaries,
    dict
  )
  chVersions = chVersions.mix(gatherPileupSummaries.out.versions)

  /*
   * CALCULATE CONTAMINATION
   */

  calculateContamination(
    bam.join(gatherPileupSummaries.out.mergedPileupFileCh)
  )
  chVersions = chVersions.mix(calculateContamination.out.versions)

  /*
   * FILTER MUTECT CALL
   */

  concatVCF.out.vcf.map{
      meta, variantCaller, vcf, index -> [meta, vcf, index]
  }.join(mergeMutect2Stats.out.mergedStatsFile)
  .set{mutect2CallsToFilter}

  mutect2CallsToFilter = params.skipMutectContamination ? 
    mutect2CallsToFilter.combine(Channel.from('NO_FILE') : 
    mutect2CallsToFilter.join(calculateContamination.out.contaminationTable)

  filterMutect2Calls(
    mutect2CallsToFilter,
    dict,
    fasta,
    fai,
    germlineResource,
    germlineResourceIndex,
    intervals
  )
  chVersions = chVersions.mix(filterMutect2Calls.out.versions)

  emit:
  versions = chVersions
  vcfUnfiltered = mutect2.out.vcf
  vcfFiltered = filterMutect2Calls.out.filteredMutect2OutputCh
  stats = mergeMutect2Stats.out.mergedStatsFile
}
