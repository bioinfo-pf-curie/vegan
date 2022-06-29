/*
 * HaplotypeCaller Flow
 */

include { haplotypeCaller } from '../../common/process/gatk/haplotypeCaller'
include { genotypeGVCFs } from '../../common/process/gatk/genotypeGVCFs'

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

  emit:
  gvcf = haplotypeCaller.out.gvcf
  vcf = genotypeGVCFs.out.vcf
  versions = chVersions
}
