/*
 * HaplotypeCaller Flow
 */

include { haplotypeCaller } from '../../common/process/gatk/haplotypeCaller'
include { mergeVCFs } from '../../common/process/gatk/mergeVCFs'

workflow haplotypeCallerSingleSampleFlow {

  take:
  bam
  intervals
  dbsnp
  dbsnpIndex
  fasta
  fastaFai
  dict

  main:
  chVersions = Channel.empty()

  // Add intervals size for groupKey function
  chIntervalsNum = intervals.collect().map{it -> [it, it.size() ]}.transpose()

  if (params.noIntervals){
    chBamIntervals = bam.map{ meta, bam, bai ->
      def newMeta = meta.clone()
      newMeta.numIntervals = 0
      [newMeta,bam,bai,[]]
    }
  }else{
    chBamIntervals = bam.combine(chIntervalsNum)
      .map { meta, bam, bai, intervals, num ->
        def newMeta = meta.clone()
        newMeta.numIntervals = num
        newMeta.intervals = intervals.getName()
        [newMeta, bam, bai, intervals]
      }
  }

  haplotypeCaller(
    chBamIntervals,
    dbsnp,
    dbsnpIndex,
    fasta,
    fastaFai,
    dict
  )
  chVersions = chVersions.mix(haplotypeCaller.out.versions)

  chVcfMerge = haplotypeCaller.out.vcf
    .join(haplotypeCaller.out.tbi)
    .map{ meta, vcf, tbi ->
      def newMeta = meta.clone()
      newMeta.remove('intervals')
      [groupKey(newMeta, newMeta.numIntervals), vcf, tbi]
    }.groupTuple()
    .branch {
      single: it[0].numIntervals <= 1
      multiple: it[0].numIntervals > 1
    }

  mergeVCFs(
    chVcfMerge.multiple,
    dict
  )
  chVersions = chVersions.mix(mergeVCFs.out.versions)

  chVcf = mergeVCFs.out.vcf
    .mix(chVcfMerge.single)
    .map{meta, vcf, tbi ->
      def newMeta = [ id: meta.id, status: meta.status, sex: meta.sex ]
      [ newMeta, vcf, tbi ]
    }

  emit:
  vcf = chVcf
  vcfFiltered = Channel.empty()
  versions = chVersions
}
