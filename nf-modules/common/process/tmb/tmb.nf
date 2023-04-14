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
  path (bed)
  val (effgsize)

  output:
  tuple val(meta), path("*_tmb.txt"), emit: tmb
  path "versions.txt"               , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def bedOpts = bed ? "--bed ${bed}" : "--effGenomeSize ${effgsize}"
  """
  pyTMB.py -i $vcf \
      ${args} \
      --dbConfig ${dbconfig} \
      --varConfig ${varconfig} \
      ${bedOpts} > ${prefix}_tmb.txt


  pyTMB.py --version > versions.txt
  """
}
