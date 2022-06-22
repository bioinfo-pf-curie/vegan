/*
 * Germine variant calling with Mutect2
 */

process mutect2 {
  tag "$meta.id"
  label 'gatk'
  label 'medMem'
  label 'lowCpu'

  input:
  tuple val(meta), path(bam), path(bai)
  path bed
  path fasta
  path fai
  path dict
  path germlineResource
  path germlineResourceIndex
  path panelOfNormals
  path panelOfNormalsIndex

  output:
  tuple val(meta), path("*.vcf.gz")     , emit: vcf
  tuple val(meta), path("*.tbi")        , emit: tbi
  tuple val(meta), path("*.stats")      , emit: stats
  path "versions.txt"                   , emit: versions

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def inputs = bam.collect{ "--input $it"}.join(" ")
  def intervalCmd = bed ? "--intervals $bed" : ""
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

  echo "GATK "\$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//' > versions.txt
  """
}
