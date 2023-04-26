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

  // Add intervals size for groupKey function
  chIntervalsNum = intervals.collect().map{it -> [it, it.size() ]}.transpose()

  // Add numIntervals in meta
  if (params.noIntervals){
    chBamIntervals = bam.map{ meta,tumorBam,tumorBai ->
      def newMeta = meta.clone()
      newMeta.numIntervals = 0
      [newMeta,tumorBam,tumorBai,[]]
    }
  }else{
    chBamIntervals = bam.combine(chIntervalsNum)
      .map { meta,tumorBam,tumorBai,intervals,num ->
        def newMeta = meta.clone()
        newMeta.numIntervals = num
        [newMeta,tumorBam,tumorBai,intervals]
      }
    }

  mutect2Tumor(
    chBamIntervals,
    fasta,
    fai,
    dict,
    germlineResource,
    germlineResourceIndex,
    panelsOfNormals,
    panelsOfNormalsIndex
  )
  chVersions = chVersions.mix(mutect2Tumor.out.versions)

  chMutect2Vcfs =  mutect2Tumor.out.vcf
    .map{ meta, vcf, tbi ->
      [groupKey(meta, meta.numIntervals), vcf, tbi]
    }.groupTuple()
    .branch {
      single: it[0].numIntervals <= 1
      multiple: it[0].numIntervals > 1
    }

  mergeVCFs(
    chMutect2Vcfs.multiple,
    dict
  )
  chVersions = chVersions.mix(mergeVCFs.out.versions)
  chMutect2MergedVcf = mergeVCFs.out.vcf.mix(chMutect2Vcfs.single)

  chMutect2Stats = mutect2Tumor.out.stats
    .map{ meta, stats ->
      [groupKey(meta, meta.numIntervals), stats]
    }.groupTuple()
    .branch {
      single: it[0].numIntervals <= 1
      multiple: it[0].numIntervals > 1
    }

  mergeMutect2Stats(
    chMutect2Stats.multiple
  )
  chVersions = chVersions.mix(mergeMutect2Stats.out.versions)
  chMutect2MergedStats = mergeMutect2Stats.out.stats.mix(chMutect2Stats.single)

  /*
   * STRAND BIAS
   */

  chMutect2LoM = mutect2Tumor.out.f1r2
    .map{ meta, fr ->
      [groupKey(meta, meta.numIntervals), fr]
    }.groupTuple()

  learnReadOrientationModel(
    chMutect2LoM
  )
  chVersions = chVersions.mix(learnReadOrientationModel.out.versions)

  /*
   * CALCULATE CONTAMINATION
   */

  getPileupSummaries(
    chBamIntervals,
    pileupSum,
    pileupSumIndex
  )
  chVersions = chVersions.mix(getPileupSummaries.out.versions)

  chPileupSum = getPileupSummaries.out.table
    .map{ meta, fr ->
      [groupKey(meta, meta.numIntervals), fr]
    }.groupTuple()
    .branch {
      single: it[0].numIntervals <= 1
      multiple: it[0].numIntervals > 1
    }

  gatherPileupSummaries(
    chPileupSum.multiple,
    dict
  )
  chVersions = chVersions.mix(gatherPileupSummaries.out.versions)
  chMergedPileupSum = gatherPileupSummaries.out.table.mix(chPileupSum.single)

  calculateContamination(
    chMergedPileupSum.map{meta,pileup->[meta, pileup, []]}
  )
  chVersions = chVersions.mix(calculateContamination.out.versions)

  /*
   * FILTER MUTECT CALL
   */

  chConta = !params.skipMutectContamination ? calculateContamination.out.contamination : chMutect2MergedVcf.map{it->[it[0],[]]}
  chContaSegment = !params.skipMutectContamination ? calculateContamination.out.segmentation : chMutect2MergedVcf.map{it->[it[0],[]]}
  chOrientationModel = !params.skipMutectOrientationModel ? learnReadOrientationModel.out.orientation :	chMutect2MergedVcf.map{it->[it[0],[]]}
  chMutect2CallsToFilter = chMutect2MergedVcf
      .join(chMutect2MergedStats)
      .join(chOrientationModel)
      .join(chConta)
      .join(chContaSegment)
      .map{meta, vcf, index, stats, orientation, conta, segment ->
        newMeta = [tumor_id:meta.tumor_id, normal_id:meta.normal_id, pair_id:meta.pair_id, id:meta.id, status:meta.status, sex:meta.sex]
        [newMeta, vcf, index, stats, orientation, conta, segment]
      }

  filterMutect2Calls(
    chMutect2CallsToFilter,
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
  vcfRaw = chMutect2MergedVcf
  vcfFiltered = bcftoolsNorm.out.vcf
  stats = chMutect2MergedStats
}
