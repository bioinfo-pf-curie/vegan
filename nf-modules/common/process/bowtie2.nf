/*
 * Alignment on reference genome with Bowtie2
 */

process bowtie2{
  tag "${meta}"
  label 'bowtie2'
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
  inputOpts = meta.singleEnd ? "-U ${reads[0]}" : "-1 ${reads[0]} -2 ${reads[1]}"
  """
  localIndex=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
  refName=`basename \${localIndex}`
  echo \$(bowtie2 --version | awk 'NR==1{print "bowtie2 "\$3}') > versions.txt
  bowtie2 -p ${task.cpus} \
          ${args} \
           -x \${localIndex} \
          $inputOpts > ${meta.id}_\${refName}.bam 2> ${meta.id}_bowtie2.log
  """
}


