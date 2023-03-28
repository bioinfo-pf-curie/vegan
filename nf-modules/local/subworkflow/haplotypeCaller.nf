/*
 * HaplotypeCaller Flow
 */

include { haplotypeCaller } from '../../common/process/gatk/haplotypeCaller'
include { genotypeGVCFs } from '../../common/process/gatk/genotypeGVCFs'
include { mergeVCFs } from '../../common/process/gatk/mergeVCFs'
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

  chGVCF = chBamIntervals
    .join(haplotypeCaller.out.gvcf)
    .map{meta,bam,bai,intervals,gvcf,index -> [meta,gvcf,index,intervals]}

  genotypeGVCFs(
    chGVCF,
    dbsnp,
    dbsnpIndex,
    fasta,
    fastaFai,
    dict
  )
  chVersions = chVersions.mix(genotypeGVCFs.out.versions)

  chVcfMerge = genotypeGVCFs.out.vcf
    .map{ meta, vcf, tbi ->
      [groupKey(meta, meta.numIntervals), vcf, tbi]
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

  bcftoolsNorm(
    chVcf,
    fasta
  )
  chVersions = chVersions.mix(bcftoolsNorm.out.versions)

  emit:
  vcf = chVcf
  vcfNorm = bcftoolsNorm.out.vcf
  versions = chVersions
}
