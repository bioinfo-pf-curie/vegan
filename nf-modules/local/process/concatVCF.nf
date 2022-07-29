/*
 * merge vcf & apply target to somatic / germline variant calling files
 */

process concatVCF {
  tag "${prefix}"
  label 'bcftools'
  label 'highCpu'
  label 'medMem'

  input:
  tuple val(meta), path(vcf)
  path(targetBed)
  path(fasta)
  path(fastaFai)

  output:
  tuple val(meta), path("*concat.vcf.gz"), path("*concat.vcf.gz.tbi"), emit: vcf
  path("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  apConcatenateVCFs.sh -g ${fasta} -i ${fastaFai} -c ${task.cpus} -o ${prefix}_concat.vcf ${args}
  bcftools --version &> versions.txt 2>&1 || true
  """
}
