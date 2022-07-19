/*
 * HaplotypeCaller Flow
 */

include { haplotypeCaller } from '../../common/process/gatk/haplotypeCaller'
include { genotypeGVCFs } from '../../common/process/gatk/genotypeGVCFs'
include { concatVCF } from '../../local/process/concatVCF'
include { collectVCFmetrics } from '../../local/process/collectVCFmetrics'

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

  emit:
  gvcf = haplotypeCaller.out.gvcf
  vcf = concatVCF.out.vcf
  mqc = collectVCFmetrics.out.mqc
  versions = chVersions
}
