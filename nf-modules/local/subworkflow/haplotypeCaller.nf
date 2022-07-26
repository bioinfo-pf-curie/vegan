/*
 * HaplotypeCaller Flow
 */

include { haplotypeCaller } from '../../common/process/gatk/haplotypeCaller'
include { genotypeGVCFs } from '../../common/process/gatk/genotypeGVCFs'
include { concatVCF } from '../../local/process/concatVCF'
include { collectVCFmetrics } from '../../local/process/collectVCFmetrics'
include { bcftoolsNorm } from '../../common/process/bcftools/bcftoolsNorm'
include { computeTransition } from '../../local/process/computeTransition'

workflow haplotypeCallerFlow {

  take:
  bqsrBam
  bed
  dbsnp
  dbsnpIndex
  fasta
  fastaFai
  dict

  main:
  chVersions = Channel.empty()

  haplotypeCaller(
    bqsrBam,
    bed,
    dbsnp,
    dbsnpIndex,
    fasta,
    fastaFai,
    dict
  )
  chVersions = chVersions.mix(haplotypeCaller.out.versions)

  genotypeGVCFs(
    haplotypeCaller.out.gvcf,
    bed,
    dbsnp,
    dbsnpIndex,
    fasta,
    fastaFai,
    dict
  )
  chVersions = chVersions.mix(genotypeGVCFs.out.versions)

  concatVCF(
    genotypeGVCFs.out.vcf,
    bed,
    fasta,
    fastaFai
  )

  chVersions = chVersions.mix(concatVCF.out.versions)

  collectVCFmetrics(
    concatVCF.out.vcf
  )

  bcftoolsNorm(
    concatVCF.out.vcf,
    fasta
    )

  bcftoolsNorm.out.vcf.view()

  computeTransition(
    bcftoolsNorm.out.vcf
    )

  emit:
  gvcf = haplotypeCaller.out.gvcf
  vcf = concatVCF.out.vcf
  vcfNorm = bcftoolsNorm.out.vcf
  transition = computeTransition.out.metrics
  mqc = collectVCFmetrics.out.mqc
  versions = chVersions
}
