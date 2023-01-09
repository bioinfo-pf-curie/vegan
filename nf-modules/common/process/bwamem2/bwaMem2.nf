/*
 * Alignement on reference genome with Bwa-mem
 */

process bwaMem2{
  tag "${meta.id}"
  label 'bwaMem2'
  label 'highCpu'
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

  bwa-mem2 \
        mem \
        -k 19 -T 30 -M \
        -t $task.cpus \
        \${localIndex} \
        $reads > ${prefix}_\${refName}.sam

  getBWAstats.sh -i ${prefix}_\${refName}.bam -p ${task.cpus} > ${prefix}_bwa.log
  echo "Bwa-mem2 "\$(bwa-mem2 version 2>&1) &> versions.txt
  """
}
