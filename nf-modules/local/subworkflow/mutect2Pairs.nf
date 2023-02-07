/*
 * Mutect2 Tumor/Normal
 */

include { mutect2 as mutect2Pairs } from '../../common/process/gatk/mutect2'
include { mergeMutect2Stats } from '../../common/process/gatk/mergeMutect2Stats'
include { learnReadOrientationModel } from '../../common/process/gatk/learnReadOrientationModel'
include { getPileupSummaries as getPileupSummariesTumor } from '../../common/process/gatk/getPileupSummaries'
include { getPileupSummaries as getPileupSummariesNormal } from '../../common/process/gatk/getPileupSummaries'
include { gatherPileupSummaries as gatherPileupSummariesTumor } from '../../common/process/gatk/gatherPileupSummaries'
include { gatherPileupSummaries as gatherPileupSummariesNormal } from '../../common/process/gatk/gatherPileupSummaries'
include { calculateContamination } from '../../common/process/gatk/calculateContamination'
include { filterMutect2Calls } from '../../common/process/gatk/filterMutect2Calls'
include { bcftoolsNorm } from '../../common/process/bcftools/bcftoolsNorm'
include { mergeVCFs } from '../../common/process/gatk/mergeVCFs'

workflow mutect2PairsFlow {

  take:
  bam // [meta, tumor_bam, tumor_bai, normal_bam, normal_bai]
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

  // Mutect2 tumor/bam inputs
  //[meta][tumor_bam, normal_bam],[tumor_bai, normal_bai]
  bam
    .map{ it -> [it[0], [it[1], it[3]], [it[2], it[4]]] }
    .set{ chBamPair }
  bam
    .map{ it -> [it[0], it[1],it[2]]}
    .set{ chBamTumor }
  bam
    .map{ it -> [it[0], it[3],it[4]]}
    .set{ chBamNormal }

  /*
   * MUTECT2 CALLS
   */

  mutect2Pairs(
    chBamPair.combine(intervals),
    fasta,
    fai,
    dict,
    germlineResource,
    germlineResourceIndex,
    panelsOfNormals,
    panelsOfNormalsIndex
  )
  chVersions = chVersions.mix(mutect2Pairs.out.versions)

  mergeVCFs(
    mutect2Pairs.out.vcf.map{it -> [it[0],it[1]]}.groupTuple(),
    dict
  )
  chVersions = chVersions.mix(mergeVCFs.out.versions)
  chMutect2Vcf = mergeVCFs.out.vcf

  mergeMutect2Stats(
    mutect2Pairs.out.stats.groupTuple()
  )
  chVersions = chVersions.mix(mergeMutect2Stats.out.versions)

  /*
   * STRAND BIAS
   */

  learnReadOrientationModel(
    mutect2Pairs.out.f1r2.groupTuple()
  )
  chVersions = chVersions.mix(learnReadOrientationModel.out.versions)

  /*
   * CALCULATE CONTAMINATION
   */

  getPileupSummariesTumor(
    chBamTumor.combine(intervals),
    pileupSum,
    pileupSumIndex
  )
  chVersions = chVersions.mix(getPileupSummariesTumor.out.versions)

  gatherPileupSummariesTumor(
    getPileupSummariesTumor.out.table.groupTuple(),
    dict
  )
  chVersions = chVersions.mix(gatherPileupSummariesTumor.out.versions)

  getPileupSummariesNormal(
    chBamNormal.combine(intervals),
    pileupSum,
    pileupSumIndex
  )
  chVersions = chVersions.mix(getPileupSummariesNormal.out.versions)

  gatherPileupSummariesNormal(
    getPileupSummariesNormal.out.table.groupTuple(),
    dict
  )
  chVersions = chVersions.mix(gatherPileupSummariesNormal.out.versions)

  calculateContamination(
    gatherPileupSummariesTumor.out.table.join(gatherPileupSummariesNormal.out.table)
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
    fai
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
  vcfFiltered = bcftoolsNorm.out.vcf
  stats = mergeMutect2Stats.out.stats
}
