/*
 * GATK collectWGSmetrics
 */

process collectWgsMetrics {
  tag "${meta.id}"
  label 'gatk'
  label 'minCpu'
  label 'lowMem'

  input:
  tuple val(meta), path(bamFiltered), path(baiFiltered)
  path(bed)
  path(fasta)
  path(dict)

  output:
  path("*metrics.txt"), emit: metrics
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def preproc = bed ? "gatk BedToIntervalList -I ${bed} -O intervals.bed -SD ${dict}" : ""
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  def args2 = task.ext.args2 ?: ''
  """
  ${preproc}
  gatk CollectWgsMetrics --help &> versions.txt 2>&1 || true
  gatk --java-options -Xmx${task.memory.toGiga()}g \
       ReorderSam \
       ${args} \
       -I ${bamFiltered} \
       -O ${bamFiltered.baseName}_reorder.bam \
       -SD ${dict} \
       --TMP_DIR ${params.gatkTmpDir} \

  gatk --java-options -Xmx${task.memory.toGiga()}g \
       CollectWgsMetrics \
       -I ${bamFiltered.baseName}_reorder.bam \
       -O ${bamFiltered.baseName}_collect_wgs_metrics.txt \
       -R ${fasta} \
       ${args2}
  """
}
