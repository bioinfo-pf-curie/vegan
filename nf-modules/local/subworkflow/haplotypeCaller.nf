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
  chBamIntervals = bam.combine(intervals)

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

  mergeVCFs(
    genotypeGVCFs.out.vcf.map{it -> [it[0],it[1]]}.groupTuple(),
    dict
  )
  chVersions = chVersions.mix(mergeVCFs.out.versions)

  bcftoolsNorm(
    mergeVCFs.out.vcf,
    fasta
  )
  chVersions = chVersions.mix(bcftoolsNorm.out.versions)

  emit:
  vcf = mergeVCFs.out.vcf
  vcfNorm = bcftoolsNorm.out.vcf
  versions = chVersions
}
