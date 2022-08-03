/*
 * HaplotypeCaller Flow
 */

include { haplotypeCaller } from '../../common/process/gatk/haplotypeCaller'
include { genotypeGVCFs } from '../../common/process/gatk/genotypeGVCFs'
//include { concatVCF } from '../../local/process/concatVCF'
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

  // concatVCF(
  //   genotypeGVCFs.out.vcf,
  //   bed,
  //   fasta,
  //   fastaFai
  // )

  //chVersions = chVersions.mix(concatVCF.out.versions)

  genotypeGVCFs.out.vcf
    .map{ it -> [it[0], it[1], it[2], [], [], []] }
    .set{ chGenoVCF }

  chGenoVCF.map{ it -> [it[0], it[1], it[2]]}.set { chGenoSimple }

  collectVCFmetrics(
    chGenoVCF
  )


  bcftoolsNorm(
    chGenoSimple,
    fasta
    )


  computeTransition(
    bcftoolsNorm.out.vcf
    )

  emit:
  gvcf = haplotypeCaller.out.gvcf
  vcf = genotypeGVCFs.out.vcf
  vcfNorm = bcftoolsNorm.out.vcf
  transition = computeTransition.out.metrics
  mqc = collectVCFmetrics.out.mqc
  versions = chVersions
}
