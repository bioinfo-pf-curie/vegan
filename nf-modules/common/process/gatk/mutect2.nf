/*
 * Somatic variant calling with Mutect2
 */

process mutect2 {
  tag "${meta.id}"
  label 'gatk'
  label 'medMem'
  label 'medCpu'

  input:
  tuple val(meta), path(bam), path(bai), path(intervals)
  path fasta
  path fai
  path dict
  path germlineResource
  path germlineResourceIndex
  path panelOfNormals
  path panelOfNormalsIndex

  output:
  tuple val(meta), path("*.vcf.gz"), path("*.tbi"), emit: vcf
  tuple val(meta), path("*.stats")                , emit: stats
  tuple val(meta), path("*f1r2.tar.gz")           , emit: f1r2, optional: true
  path("versions.txt")                            , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def inputs = bam.collect{ "--input $it"}.join(" ")
  def intervalCmd = intervals ? "--intervals $intervals" : ""
  def ponCmd = panelOfNormals ? "--panel-of-normals $panelOfNormals" : ""
  def grCmd = germlineResource ? "--germline-resource $germlineResource" : ""

  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" Mutect2 \\
       $inputs \\
       --output ${prefix}.vcf.gz \\
       --reference $fasta \\
       $ponCmd \\
       $grCmd \\
       $intervalCmd \\
       --tmp-dir . \\
       $args

  echo "GATK "\$(gatk --version 2>&1 | grep \\(GATK\\) | sed 's/^.*(GATK) v//') > versions.txt
  """
}
