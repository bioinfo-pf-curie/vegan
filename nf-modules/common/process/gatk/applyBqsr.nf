/*
 * Apply BQSR - detects systematic errors made by the sequencing machine
 */

process applyBQSR {
  tag "${meta.id}"
  label 'gatk'
  label 'lowMem'
  label 'lowCpu'

  input:
  tuple val(meta), path(bamFiltered), path(bamFilteredBai), path(recalTable), path(intervals)
  path(fasta)
  path(fastaFai)
  path(dict)

  output:
  tuple val(meta),path("*.recal.bam"), emit:bam
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def intervalsCmd = intervals ? "--intervals $intervals" : ""
  """
  gatk --java-options -Xmx${task.memory.toGiga()}g \
       ApplyBQSR \
       -R ${fasta} \
       --input ${bamFiltered} \
       --output ${prefix}.recal.bam \
       ${intervalsCmd} \
       ${args} \
       --bqsr-recal-file ${recalTable}
  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
 }
