/*
 * GATK MarkDuplicates
 */

process markDuplicates {
  tag "${meta.id}"
  label 'gatk'
  label 'minCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bam), path(bai)
  path(fasta)
  path(fai)

  output:
  tuple val(meta), path('*markDups.bam'), optional:true, emit: bam
  tuple val(meta), path('*markDups.cram'), optional:true, emit: cram
  tuple val(meta), path("*.metrics"), emit: metrics
  path('versions.txt'), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  def reference = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""
  def extension = bam.getExtension()
  def availMem = (task.memory.mega*0.8).intValue()
  """
  gatk --java-options "-Xmx${availMem}M -XX:-UsePerfData" \\
    MarkDuplicates \\
    --INPUT ${bam} \\
    --OUTPUT ${prefix}_markDups.bam \\
    --METRICS_FILE ${prefix}.metrics \\
    --TMP_DIR . \\
    ${reference} \\
    $args

  # If cram files are wished as output, the run samtools for conversion
  if [[ ${extension} == 'cram' ]]; then
      samtools view -Ch -T ${fasta} -o ${prefix}_markDups.cram ${prefix}_markDups.bam
      rm ${prefix}_markDups.bam
  fi
 
  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
