/*
 * HaplotypeCaller Flow
 */

include { haplotypeCaller } from '../../local/process/haplotypeCaller'
include { genotypeGVCFs } from '../../local/process/genotypeGVCFs'

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
      bed.collect(),
      dbsnp.collect(),
      dbsnpIndex.collect(),
      fasta.collect(),
      fastaFai.collect(),
      dict.collect()
      )

    genotypeGVCFs(
      haplotypeCaller.out.hcGvcf,
      bed.collect(),
      dbsnp.collect(),
      dbsnpIndex.collect(),
      fasta.collect(),
      fastaFai.collect(),
      dict.collect()
      )

  emit:
    hcGvcf = haplotypeCaller.out.hcGvcf
    hcVcf = genotypeGVCFs.out.hcVcf

}
