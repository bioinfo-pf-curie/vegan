/*
 * Mutect2 Tumor/Normal
 */

include { mutect2 as mutect2Pairs } from '../../common/process/gatk/mutect2'
include { mergeMutect2Stats } from '../../common/process/gatk/mergeMutect2Stats'
//include { concatVCF } from '../../local/process/concatVCF'
include { learnReadOrientationModel } from '../../common/process/gatk/learnReadOrientationModel'
include { getPileupSummaries as getPileupSummariesTumor } from '../../common/process/gatk/getPileupSummaries'
include { getPileupSummaries as getPileupSummariesNormal } from '../../common/process/gatk/getPileupSummaries'
include { gatherPileupSummaries as gatherPileupSummariesTumor } from '../../common/process/gatk/gatherPileupSummaries'
include { gatherPileupSummaries as gatherPileupSummariesNormal } from '../../common/process/gatk/gatherPileupSummaries'
include { calculateContamination } from '../../common/process/gatk/calculateContamination'
include { filterMutect2Calls } from '../../common/process/gatk/filterMutect2Calls'
include { collectVCFmetrics } from '../../local/process/collectVCFmetrics'
include { bcftoolsNorm } from '../../common/process/bcftools/bcftoolsNorm'
include { computeTransition } from '../../local/process/computeTransition'

workflow mutect2PairsFlow {

  take:
  bam
  bed
  fasta
  fai
  dict
  germlineResource
  germlineResourceIndex
  pileupSum
  pileupSumIndex
  panelsOfNormals
  panelsOfNormalsIndex
  intervals

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
    chBamPair,
    bed,
    fasta,
    fai,
    dict,
    germlineResource,
    germlineResourceIndex,
    panelsOfNormals,
    panelsOfNormalsIndex
  )
  chVersions = chVersions.mix(mutect2Pairs.out.versions)

  mergeMutect2Stats(
    mutect2Pairs.out.stats
  )

  // concatVCF(
  //   mutect2.out.vcf,
  //   bed,
  //   fasta,
  //   fai
  // )
  // chVersions = chVersions.mix(concatVCF.out.versions)

  /*
   * STRAND BIAS
   */

   learnReadOrientationModel(
     mutect2Pairs.out.f1r2
   )

  /*
   * CALCULATE CONTAMINATION
   */

  getPileupSummariesTumor(
    chBamTumor,
    intervals,
    pileupSum,
    pileupSumIndex,
    bed
  )
  chVersions = chVersions.mix(getPileupSummariesTumor.out.versions)

  gatherPileupSummariesTumor(
    getPileupSummariesTumor.out.table,
    dict
  )
  chVersions = chVersions.mix(gatherPileupSummariesTumor.out.versions)

  getPileupSummariesNormal(
    chBamNormal,
    intervals,
    pileupSum,
    pileupSumIndex,
    bed
  )
  chVersions = chVersions.mix(getPileupSummariesNormal.out.versions)

  gatherPileupSummariesNormal(
    getPileupSummariesNormal.out.table,
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
    mutect2Pairs.out.vcf
      .join(mergeMutect2Stats.out.stats)
      .join(learnReadOrientationModel.out.orientation)
      .map{meta, vcf, index, stats, orientation -> [meta, vcf, index, stats, orientation, [], []] }
      .set { mutect2CallsToFilter }
  }else{
    mutect2Pairs.out.vcf
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

  /*
   * VCF normalisation
   */

  collectVCFmetrics(
    filterMutect2Calls.out.vcf,
  )

  filterMutect2Calls.out.vcf
  .map{ it -> [it[0], it[3], it[4]] }
  .set{ chFiltSimple }

  bcftoolsNorm(
    chFiltSimple,
    fasta
  )

  computeTransition(
    bcftoolsNorm.out.vcf
  )
  chVersions = chVersions.mix(bcftoolsNorm.out.versions)

  emit:
  versions = chVersions
  vcfUnfiltered = mutect2Pairs.out.vcf
  vcfFiltered = filterMutect2Calls.out.vcf
  vcfFilteredNorm = bcftoolsNorm.out.vcf
  transition = computeTransition.out.metrics
  mqc = collectVCFmetrics.out.mqc
  stats = mergeMutect2Stats.out.stats
}
