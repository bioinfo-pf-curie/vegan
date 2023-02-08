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

  chBamPairIntervals = params.noIntervals ? chBamPair.map{meta,bam,bai -> [meta,bam,bai,[]]} : chBamPair.combine(intervals)

  mutect2Pairs(
    chBamPairIntervals,
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
  chMutect2Vcf = params.noIntervals || params.targetBed ? mutect2Pairs.out.vcf : mergeVCFs.out.vcf

  mergeMutect2Stats(
    mutect2Pairs.out.stats.groupTuple()
  )
  chVersions = chVersions.mix(mergeMutect2Stats.out.versions)
  chMutect2Stats = params.noIntervals || params.targetBed ? mutect2Pairs.out.stats : mergeMutect2Stats.out.stats

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

  chBamTumorIntervals = params.noIntervals ? chBamTumor.map{meta,bam,bai -> [meta,bam,bai,[]]} : chBamTumor.combine(intervals)

  getPileupSummariesTumor(
    chBamTumorIntervals,
    pileupSum,
    pileupSumIndex
  )
  chVersions = chVersions.mix(getPileupSummariesTumor.out.versions)

  gatherPileupSummariesTumor(
    getPileupSummariesTumor.out.table.groupTuple(),
    dict
  )
  chVersions = chVersions.mix(gatherPileupSummariesTumor.out.versions)
  chPileupSumTumor = params.noIntervals || params.targetBed ? getPileupSummariesTumor.out.table : gatherPileupSummariesTumor.out.table

  chBamNormalIntervals = params.noIntervals ? chBamNormal.map{meta,bam,bai -> [meta,bam,bai,[]]} : chBamNormal.combine(intervals)

  getPileupSummariesNormal(
    chBamNormalIntervals,
    pileupSum,
    pileupSumIndex
  )
  chVersions = chVersions.mix(getPileupSummariesNormal.out.versions)

  gatherPileupSummariesNormal(
    getPileupSummariesNormal.out.table.groupTuple(),
    dict
  )
  chVersions = chVersions.mix(gatherPileupSummariesNormal.out.versions)
  chPileupSumNormal = params.noIntervals || params.targetBed ? getPileupSummariesNormal.out.table : gatherPileupSummariesNormal.out.table

  calculateContamination(
    chPileupSumTumor.join(chPileupSumNormal)
  )
  chVersions = chVersions.mix(calculateContamination.out.versions)

  /*
   * FILTER MUTECT CALL
   */

  if (params.skipMutectContamination){
    chMutect2Vcf
      .join(chMutect2Stats)
      .join(learnReadOrientationModel.out.orientation)
      .map{meta, vcf, index, stats, orientation -> [meta, vcf, index, stats, orientation, [], []] }
      .set { mutect2CallsToFilter }
  }else{
    chMutect2Vcf
      .join(chMutect2Stats)
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
  stats = chMutect2Stats
}
