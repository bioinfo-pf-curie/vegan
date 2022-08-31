/*
 * Somatic variant calling with Mutect2
 */

process mutect2 {
  tag "${meta.tumor_id}_vs_${meta.normal_id}"
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
  tuple val(meta), path("*.vcf.gz"), path("*.tbi"), emit: vcf
  tuple val(meta), path("*.stats")                , emit: stats
  tuple val(meta), path("*f1r2.tar.gz")           , emit: f1r2, optional: true
  path("versions.txt")                            , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.tumor_id}_vs_${meta.normal_id}"
  def inputs = bam.collect{ "--input $it"}.join(" ")
  def intervalCmd = bed ? "--intervals $bed" : ""
  def ponCmd = panelOfNormals ? "--panel-of-normals $panelOfNormals" : ""
  def grCmd = germlineResource ? "--germline-resource $germlineResource" : ""

  """
  gatk --java-options "-Xmx${task.memory.toGiga()}g" Mutect2 \\
       $inputs \\
       --output ${prefix}_Mutect2_unfiltered.vcf.gz \\
       --reference $fasta \\
       $ponCmd \\
       $grCmd \\
       $intervalCmd \\
       --tmp-dir . \\
       $args

  echo "GATK "\$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//' > versions.txt

  """
}
