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
include { collectVCFmetrics } from '../../local/process/collectVCFmetrics'
include { bcftoolsNorm } from '../../common/process/bcftools/bcftoolsNorm'
include { computeTransition } from '../../local/process/computeTransition'

workflow mutect2TumorOnlyFlow {

  take:
  bam //[meta][tumor_bam, tumor_bai]
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

  /*
   * MUTECT2 CALLS
   */

  mutect2Tumor(
    bam,
    bed,
    fasta,
    fai,
    dict,
    germlineResource,
    germlineResourceIndex,
    panelsOfNormals,
    panelsOfNormalsIndex
  )
  chVersions = chVersions.mix(mutect2Tumor.out.versions)

  mergeMutect2Stats(
    mutect2Tumor.out.stats
  )

  learnReadOrientationModel(
     mutect2Tumor.out.f1r2
  )

  /*
   * CALCULATE CONTAMINATION
   */

  getPileupSummaries(
    bam,
    intervals,
    pileupSum,
    pileupSumIndex,
    bed
  )
  chVersions = chVersions.mix(getPileupSummaries.out.versions)

  gatherPileupSummaries(
    getPileupSummaries.out.table,
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
    mutect2Tumor.out.vcf
      .join(mergeMutect2Stats.out.stats)
      .join(learnReadOrientationModel.out.orientation)
      .map{meta, vcf, index, stats, orientation -> [meta, vcf, index, stats, orientation, [], []] }
      .set { mutect2CallsToFilter }
  }else{
    mutect2Tumor.out.vcf
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
  vcfUnfiltered = mutect2Tumor.out.vcf
  vcfFiltered = filterMutect2Calls.out.vcf
  vcfFilteredNorm = bcftoolsNorm.out.vcf
  transition = computeTransition.out.metrics
  mqc = collectVCFmetrics.out.mqc
  stats = mergeMutect2Stats.out.stats
}
