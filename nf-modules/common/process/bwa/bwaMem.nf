/*
 * Alignement on reference genome with Bwa-mem
 */

process bwaMem{
  tag "${meta.id}"
  label 'bwa'
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
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  localIndex=`find -L ./ -name "*.amb" | sed 's/.amb//'`
  refName=`basename \${localIndex}`

  bwa \
    mem \
    $args \
    -t $task.cpus \
    \${localIndex} \
    $reads > ${prefix}_\${refName}.sam 

  echo "Bwa-mem "\$(bwa 2>&1 | grep Version | cut -d" " -f2) &> versions.txt

  """
}
