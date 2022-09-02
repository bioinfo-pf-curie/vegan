/*
 * Mostdepth - calculate sequencing depth
 */

process mosdepth {
  tag "${meta.id}"
  label 'mosdepth'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bamFiltered), path(baiFiltered)
  path(bed)

  output:
  path("*.*.txt"), emit: metrics
  path("*{.bed.gz,.bed.gz.csi}"), emit: bedcov
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  bedOpts = bed ? "--by ${bed}" : ''
  """
  mosdepth --version &> versions.txt 2>&1 || true
  mosdepth -t ${task.cpus} ${args} ${bedOpts} ${prefix} ${bamFiltered}
  """
}
