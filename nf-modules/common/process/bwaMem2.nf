/*
 * Alignement on reference genome with Bwa-mem
 */

process bwaMem2{
  tag "${meta.id}"
  label 'bwaMem2'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(meta), path(reads)
  path(index)

  output:
  tuple val(meta), path("*.bam"), emit: bam
  path("*.log"), emit: logs
  path("versions.txt"), emit: versions

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
        \$INDEX \
        $reads \
        | samtools view -bS -@ $task.cpus -o ${prefix}_\${refName}.bam -

  getBWAstats.sh -i ${prefix}_\${refName}.bam -p ${task.cpus} > ${prefix}_bwa.log
  echo "Bwa-mem2 "\$(bwa-mem2 version 2>&1) &> versions.txt
  """
}
