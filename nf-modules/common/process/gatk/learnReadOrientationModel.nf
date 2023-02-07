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
  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" LearnReadOrientationModel \\
    ${inputList} \\
    --tmp-dir . \\
    --output ${prefix}_read-orientation-model.tar.gz

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
