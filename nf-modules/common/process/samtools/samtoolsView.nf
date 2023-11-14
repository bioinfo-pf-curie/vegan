/*
 * Samtools - View
 */

process samtoolsView {
  tag "${meta.id}"
  label 'samtools'
  label 'medCpu'
  label 'medMem'

  input:
  tuple val(meta), path(input)
  path(fasta)

  output:
  tuple val(meta), path ("*.sam"), optional:true, emit: sam
  tuple val(meta), path ("*.bam"), optional:true, emit: bam
  tuple val(meta), path ("*.cram"), optional:true, emit: cram
  path("versions.txt") , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def extension = args.contains("--output-fmt sam") ? "sam" :
                  args.contains("--output-fmt bam") ? "bam" :
                  args.contains("--output-fmt cram") ? "cram" :
                  input.getExtension()
  def reference = fasta ? "--reference ${fasta}" : ""
  """
  samtools view \\
    ${args} \\
    -@  ${task.cpus} \\
    -o ${prefix}.${extension} \\
    ${reference} \\
    ${input}

  echo \$(samtools --version | head -1 ) > versions.txt
  """
}
