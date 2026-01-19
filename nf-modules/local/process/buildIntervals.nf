process buildIntervals {
  label 'unix'
  label 'minCpu'
  label 'lowMem'

  input:
  path(fai)

  output:
  path("*.bed"), emit: bed

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''

  script:
  """
  awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fai} > ${fai.baseName}.bed
  """
}
