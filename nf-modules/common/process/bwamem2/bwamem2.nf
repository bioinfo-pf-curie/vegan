/*
 * Alignement on reference genome with Bwa-mem
 */

process bwamem2{
  tag "${meta.id}"
  label 'bwamem2'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(meta), path(reads)
  path(index)

  output:
  tuple val(meta), path("*.bam"), emit: bam
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  localIndex=`find -L ./ -name "*.amb" | sed 's/.amb//'`
  refName=\$(basename \${localIndex})

   bwa-mem2 \
        mem \
        $args \
        -t $task.cpus \
        \$localIndex \
        $reads \
        | samtools view -bS -@ $task.cpus -o ${prefix}_\${refName}.bam -

  echo "Bwa-mem2 "\$(bwa-mem2 version 2> /dev/null) &> versions.txt
  """
}
