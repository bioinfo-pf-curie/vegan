/*
 * GATK: Learn Read Orientation Model: to identify strand biais artefacts from FFPE or suspicious samples
 */

process learnReadOrientationModel {
  label 'gatk'
  label 'minCpu'
  label 'medMem'

  tag "${meta.tumor_id}_vs_${meta.normal_id}"

  input:
  tuple val(meta), path(f1r2)

  output:
  path("*_read-orientation-model.tar.gz"), emit: readOrientation

  when:
  task.ext.when == null || task.ext.when

  script:
  //def args = task.ext.args ?: ''
  """

  gatk --java-options "-Xmx${task.memory.toGiga()}g" \
  LearnReadOrientationModel \
    -I ${f1r2} \
    -O ${meta.tumor_id}_vs_${meta.normal_id}_read-orientation-model.tar.gz
  """
}
