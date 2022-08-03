/*
 * Tumor mutational burden
 */

process tmb {
  tag "${meta.id}"
  label 'minCpu'
  label 'lowMem'
  label 'tmb'

  input:
  tuple val(meta), path(vcf), path (dbconfig), path (varconfig)
  path (target_bed)

  output:
  tuple val(meta), path("*_tmb.txt"), emit: tmb
  path "versions.txt"           , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def target_bed = target_bed ? "--bed ${target_bed}" : ""
  """
  pyTMB.py -i $vcf \
      --sample ${meta.tumor_id} \
      --dbConfig ${dbconfig} \
      --varConfig ${varconfig} \
      ${target_bed} \
      $args > ${meta.tumor_id}"_vs"${meta.normal_id}"_tmb.txt"

  pyTMB.py --version > versions.txt
  """
}
