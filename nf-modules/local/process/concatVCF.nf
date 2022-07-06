/*
 * merge vcf & apply target to somatic / germline variant calling files
 */

process concatVCF {
  tag "${variantCaller}-${meta.id}"
  label 'bcftools'
  label 'highCpu'
  label 'medMem'

  input:
  tuple val(meta), val(variantCaller), path(vcFiles)
  path(targetBed)
  path(fasta)
  path(fastaFai)

  output:
  tuple val(meta), val(variantCaller), path("*_*.vcf.gz"), path("*_*.vcf.gz.tbi"), emit: vcf
  path("versions.txt")                   , emit: versions

  script:
  if (variantCaller == 'HaplotypeCaller')
    outputFile = "${meta.id}_HaplotypeCaller.vcf"
  else if (variantCaller == 'Mutect2')
    outputFile = "${meta.tumor_id}_vs_${meta.normal_id}_Mutect2_unfiltered.vcf"
  else
    outputFile = "${meta.id}_${variantCaller}.vcf"

  def args = task.ext.args ?: ''
  """
  apConcatenateVCFs.sh -g ${fasta} -i ${fastaFai} -c ${task.cpus} -o ${outputFile} ${args}
  bcftools --version &> versions.txt 2>&1 || true
  """
}
