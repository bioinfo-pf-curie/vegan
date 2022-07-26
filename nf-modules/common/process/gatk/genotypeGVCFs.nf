/*
 * Germine variant calling
 */

process genotypeGVCFs {
  tag "${fileID}"
  label 'gatk'

  input:
  tuple val(meta), path(gvcf)
  path(bed)
  path(dbsnp)
  path(dbsnpIndex)
  path(fasta)
  path(fastaFai)
  path(dict)

  output:
  tuple val(meta), path("${fileID}.vcf"), emit: vcf
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def args2 = task.ext.args2 ?: ''
  //prefix = task.ext.prefix ?: "${meta.id}"
  fileID = "${meta.status}" == "pair" ? "${meta.tumor_id}_vs_${meta.normal_id}" : "${meta.status}" == "tumor" ? "${meta.tumor_id}" : "${meta.normal_id}"
  """
  gatk --java-options -Xmx${task.memory.toGiga()}g \
    IndexFeatureFile -I ${gvcf}

  gatk --java-options -Xmx${task.memory.toGiga()}g \
    GenotypeGVCFs \
    -R ${fasta} \
    ${args} \
    ${args2} \
    -V ${gvcf} \
    -O ${fileID}.vcf

  echo "GATK "\$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//' > versions.txt
  """
}
