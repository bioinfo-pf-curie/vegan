/*
 * Concat multiple VCFs from the same sample
 */

process bcftoolsConcat {
  label 'bcftools'
  label 'medCpu'
  label 'medMem'
  tag "${meta.id}"

  input:
  tuple val(meta), path(vcfs), path(idx)

  output:
  tuple val(meta), path("*.vcf{.gz,.gz.tbi}"), emit: vcf
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  
  """
  bcftools concat ${args} ${vcfs} \
    --threads ${task.cpus} \
    -o ${prefix}_concat.vcf.gz
  tabix ${prefix}_concat.vcf.gz

  echo \$(bcftools --version | head -1) > versions.txt
  echo "tabix "\$(tabix 2>&1 | awk '\$1~"Version"{print \$2}') >> versions.txt
  """
}