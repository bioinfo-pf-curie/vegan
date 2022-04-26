/*
 * getSeqDepth:
 * External parameters :
 */

process getSeqDepth {
  tag "${meta.id}"
  label 'mosdepth'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bamFiltered), path(baiFiltered)
  path(bed)

  output:
  path("*.*.txt"), emit: metrics
  path("*{.bed.gz,.bed.gz.csi}"), emit: mosdepthBed
  path("versions.txt"), emit: versions

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  """
  mosdepth --version &> versions.txt 2>&1 || true
  mosdepth -t ${task.cpus} -n --quantize 0:1:10:50:100: ${args} ${prefix} ${bamFiltered}
  """
}
