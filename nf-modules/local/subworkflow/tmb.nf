/*
 * Tumor Mutational Burden
 */

include { tmb } from '../../common/process/tmb/tmb'

workflow tmbFlow {

  take:
  vcf
  bed

  main:
  chVersions = Channel.empty()

  dbConfig = file("$projectDir/assets/tmb/snpeff.yml", checkIfExists: true)
  varConfig = file("$projectDir/assets/tmb/mutect2.yml", checkIfExists: true)

  tmb(
    vcf.map{ it -> [it[0], it[1][0], dbConfig, varConfig] },
    bed
  )
  chVersions = chVersions.mix(tmb.out.versions)

  emit:
  logs = tmb.out.logs
  versions = chVersions
}
