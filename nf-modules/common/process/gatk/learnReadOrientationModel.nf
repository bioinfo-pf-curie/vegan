/*
 * GATK: Learn Read Orientation Model: to identify strand biais artefacts from FFPE or suspicious samples
 */

process learnReadOrientationModel {
  label 'gatk'
  label 'minCpu'
  label 'medMem'
  tag "${meta.id}"

  input:
  tuple val(meta), path(f1r2)

  output:
  tuple val(meta), path("*_read-orientation-model.tar.gz"), emit: orientation
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def inputList = f1r2.collect{"--input $it"}.join(' ')
  def args = task.ext.args ?: ''
  def args2 = task.ext.args2 ?: ''

  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g ${args2}" LearnReadOrientationModel \\
    ${inputList} \\
    --tmp-dir . \\
    --output ${prefix}_read-orientation-model.tar.gz \\
    ${args}

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
