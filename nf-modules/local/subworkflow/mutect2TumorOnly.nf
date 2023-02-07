/*
 * Mutect2 Tumor Only
 */

include { mutect2 as mutect2Tumor} from '../../common/process/gatk/mutect2'
include { mergeMutect2Stats } from '../../common/process/gatk/mergeMutect2Stats'
include { learnReadOrientationModel } from '../../common/process/gatk/learnReadOrientationModel'
include { getPileupSummaries } from '../../common/process/gatk/getPileupSummaries'
include { gatherPileupSummaries } from '../../common/process/gatk/gatherPileupSummaries'
include { calculateContamination } from '../../common/process/gatk/calculateContamination'
include { filterMutect2Calls } from '../../common/process/gatk/filterMutect2Calls'
include { bcftoolsNorm } from '../../common/process/bcftools/bcftoolsNorm'
include { mergeVCFs } from '../../common/process/gatk/mergeVCFs'

workflow mutect2TumorOnlyFlow {

  take:
  bam //[meta][tumor_bam, tumor_bai]
  intervals
  fasta
  fai
  dict
  germlineResource
  germlineResourceIndex
  pileupSum
  pileupSumIndex
  panelsOfNormals
  panelsOfNormalsIndex

  main:
  chVersions = Channel.empty()

  /*
   * MUTECT2 CALLS
   */

  mutect2Tumor(
    bam.combine(intervals),
    fasta,
    fai,
    dict,
    germlineResource,
    germlineResourceIndex,
    panelsOfNormals,
    panelsOfNormalsIndex
  )
  chVersions = chVersions.mix(mutect2Tumor.out.versions)

  mergeVCFs(
    mutect2Tumor.out.vcf.map{it -> [it[0],it[1]]}.groupTuple(),
    dict
  )
  chVersions = chVersions.mix(mergeVCFs.out.versions)
  chMutect2Vcf = mergeVCFs.out.vcf

  mergeMutect2Stats(
    mutect2Tumor.out.stats.groupTuple()
  )
  chVersions = chVersions.mix(mergeMutect2Stats.out.versions)

  /*
   * STRAND BIAS
   */

  learnReadOrientationModel(
     mutect2Tumor.out.f1r2
  )
  chVersions = chVersions.mix(learnReadOrientationModel.out.versions)

  /*
   * CALCULATE CONTAMINATION
   */

  getPileupSummaries(
    bam.combine(intervals),
    pileupSum,
    pileupSumIndex
  )
  chVersions = chVersions.mix(getPileupSummaries.out.versions)

  gatherPileupSummaries(
    getPileupSummaries.out.table.groupTuple(),
    dict
  )
  chVersions = chVersions.mix(gatherPileupSummaries.out.versions)

  calculateContamination(
    gatherPileupSummaries.out.table.map{meta, table -> [meta, table, []]}
  )
  chVersions = chVersions.mix(calculateContamination.out.versions)

  /*
   * FILTER MUTECT CALL
   */

  if (params.skipMutectContamination){
    chMutect2Vcf
      .join(mergeMutect2Stats.out.stats)
      .join(learnReadOrientationModel.out.orientation)
      .map{meta, vcf, index, stats, orientation -> [meta, vcf, index, stats, orientation, [], []] }
      .set { mutect2CallsToFilter }
  }else{
    chMutect2Vcf
      .join(mergeMutect2Stats.out.stats)
      .join(learnReadOrientationModel.out.orientation)
      .join(calculateContamination.out.contamination)
      .join(calculateContamination.out.segmentation)
      .set{ mutect2CallsToFilter }
  }

  filterMutect2Calls(
    mutect2CallsToFilter,
    dict,
    fasta,
    fai,
  )
  chVersions = chVersions.mix(filterMutect2Calls.out.versions)

  bcftoolsNorm(
    filterMutect2Calls.out.vcf,
    fasta
  )
  chVersions = chVersions.mix(bcftoolsNorm.out.versions)

  emit:
  versions = chVersions
  vcfRaw = chMutect2Vcf
  vcfFiltered = filterMutect2Calls.out.vcf
  stats = mergeMutect2Stats.out.stats
}
