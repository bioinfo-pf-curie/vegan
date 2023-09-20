/*
 * Alignement on reference genome with DragMap
 */

process dragmap{
  tag "${meta.id}"
  label 'dragmap'
  label 'highCpu'
  label 'extraMem'

  input:
  tuple val(meta), path(reads)
  path(hashmap)
  val(sortBam)

  output:
  tuple val(meta), path("*.bam"), emit: bam
  path("*.log"), emit: logs
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def readsCmd = meta.single_end ? "-1 $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def samtoolsCmd = sort_bam ? 'sort' : 'view'
  """
  dragen-os \\
    -r $hashmap \\
    $args \\
    --num-threads $task.cpus \\
    $readsCmd \\
    2> ${prefix}.dragmap.log \\
    | samtools ${samtoolsCmd} -bS --threads ${task.cpus} -o ${prefix}.bam -

  echo "DragMap "\$(dragen-os --version 2>&1) > versions.txt
  """
}
