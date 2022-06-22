/*
 * BQSR - detects systematic errors made by the sequencing machine
 */

process baseRecalibrator {
  tag "${meta.id}"
  label 'gatk'
  label 'minCpu'
  label 'medMem'

  input:
  tuple val(meta), path(bamFiltered), path(bamFilteredBai)
  path(bed)
  path(dbsnp)
  path(dbsnpIndex)
  path(fasta)
  path(fastaFai)
  path(knownIndels)
  path(knownIndelsIndex)
  path(dict)

  output:
  tuple val(meta), path(bamFiltered), path(bamFilteredBai), path("*.recal.table"), emit: table
  path("versions.txt"), emit: versions

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  dbsnpOptions = dbsnp.collect{"--known-sites ${it}"}.join(' ')
  knownOptions = knownIndels.collect{"--known-sites ${it}"}.join(' ')
  """
  gatk --java-options -Xmx${task.memory.toGiga()}g \
      BaseRecalibrator \
      -I ${bamFiltered} \
      -O ${prefix}.recal.table \
      --tmp-dir ${params.gatkTmpDir} \
      -R ${fasta} \
      ${args} \
      ${dbsnpOptions} \
      ${knownOptions} \
      --verbosity INFO

  echo "GATK "\$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//' > versions.txt
  """
}
