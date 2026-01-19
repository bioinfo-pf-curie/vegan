/*
 * Apply BQSR - detects systematic errors made by the sequencing machine
 */

process applyBQSR {
  tag "${meta.id}"
  label 'gatk'
  label 'lowMem'
  label 'lowCpu'

  input:
  tuple val(meta), path(bam), path(bai), path(recalTable), path(intervals)
  path(fasta)
  path(fastaFai)
  path(dict)

  output:
  tuple val(meta), path("*.recal.bam"), emit:bam, optional: true
  tuple val(meta), path("*.recal.cram"), emit: cram, optional: true
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def intervalsCmd = intervals ? "--intervals $intervals" : ""
  def extension = bam.getExtension()
  def availMem = (task.memory.mega*0.8).intValue()
  """
  gatk --java-options -Xmx${availMem}M \
       ApplyBQSR \
       --reference ${fasta} \
       --input ${bam} \
       --output ${prefix}.recal.${extension} \
       --bqsr-recal-file ${recalTable} \
       ${intervalsCmd} \
       --tmp-dir . \
       ${args}

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
 }
