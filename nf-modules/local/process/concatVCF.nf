/*
 * merge vcf & apply target to somatic / germline variant calling files
 */

process concatVCF {
  tag "${fileID}"
  label 'bcftools'
  label 'highCpu'
  label 'medMem'

  input:
  tuple val(meta), path(vcf)
  path(targetBed)
  path(fasta)
  path(fastaFai)

  output:
  tuple val(meta), path("*concat.vcf.gz*"), emit: vcf
  path("versions.txt"), emit: versions

  script:
  //def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  fileID = "${meta.status}" == "pair" ? "${meta.tumor_id}_vs_${meta.normal_id}" : "${meta.status}" == "tumor" ? "${meta.tumor_id}" : "${meta.normal_id}"
  """
  apConcatenateVCFs.sh -g ${fasta} -i ${fastaFai} -c ${task.cpus} -o ${fileID}_concat.vcf ${args}
  bcftools --version &> versions.txt 2>&1 || true
  """
}
