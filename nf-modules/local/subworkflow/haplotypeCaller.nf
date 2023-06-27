/*
 * HaplotypeCaller Flow
 */

include { haplotypeCaller } from '../../common/process/gatk/haplotypeCaller'
include { genotypeGVCFs } from '../../common/process/gatk/genotypeGVCFs'
include { mergeVCFs } from '../../common/process/gatk/mergeVCFs'
include { filterVcf } from '../../local/process/filterVcf'
include { bcftoolsNorm } from '../../common/process/bcftools/bcftoolsNorm'

workflow haplotypeCallerFlow {

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

  //chBamIntervals = params.noIntervals ? bam.map{ meta,bam,bai -> [meta,bam,bai,[]]} : bam.combine(intervals)
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

//  chHCvcf = chBamIntervals
//    .join(haplotypeCaller.out.vcf)
//    .join(haplotypeCaller.out.vcf.tbi)
//    .map{meta,bam,bai,intervals,vcf,tbi -> [meta,vcf,tbi,intervals]}.view()

//  genotypeGVCFs(
//    chGVCF,
//    dbsnp,
//    dbsnpIndex,
//    fasta,
//    fastaFai,
//    dict
//  )
//  chVersions = chVersions.mix(genotypeGVCFs.out.versions)

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

  /*
   * POST-PROCESSING
   */

  //filterVcf(
  //  chVcf
  //)

  bcftoolsNorm(
    chVcf,
    fasta
  )
  chVersions = chVersions.mix(bcftoolsNorm.out.versions)

  emit:
  vcf = chVcf
  vcfFiltered = Channel.empty() //bcftoolsNorm.out.vcf
  versions = chVersions
}
