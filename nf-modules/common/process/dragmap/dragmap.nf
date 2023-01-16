/*
 * Alignement on reference genome with DragMap
 */

process dragmap{
  tag "${meta.id}"
  label 'dragmap'
  label 'extraCpu'
  label 'extraMem'

  input:
  tuple val(meta), path(reads)
  path(index)

  output:
  tuple val(meta), path("*.sam"), emit: sam
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def readsCmd = meta.single_end ? "-1 $reads" : "-1 ${reads[0]} -2 ${reads[1]}"
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  localIndex=`find -L ./ -name "*.amb" | sed 's/.amb//'`
  refName=`basename \${localIndex}`

dragen-os \\
    -r $index \\
    $args \\
    --num-threads $task.cpus \\
    $readsCmd \\
    > ${prefix}_\${refName}.sam

  echo "DragMap "\$(dragen-os --version 2>&1) > versions.txt
  """
}
