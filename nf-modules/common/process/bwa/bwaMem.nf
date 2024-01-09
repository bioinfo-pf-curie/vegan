/*
 * Alignement on reference genome with Bwa-mem
 */

process bwaMem{
  tag "${meta.id}"
  label 'bwa'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(meta), path(reads)
  path(index)
  val(sortBam)

  output:
  tuple val(meta), path("*.bam"), emit: bam
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def samtoolsCmd = sortBam ? 'sort' : 'view'
  def bwaRGL0 = '-R "@RG\\tID:'
  def bwaRGL = "${bwaRGL0}${prefix}\\tSM:"
  def bwaRGL2 = "${bwaRGL}${prefix}"
  def bwaRGR = '\\tPL:illumina" '
  def bwaRG = "${bwaRGL2}${bwaRGR}"

  """
  localIndex=`find -L ./ -name "*.amb" | sed 's/.amb//'`
  refName=`basename \${localIndex}`

  bwa \
    mem \
    $params.bwaOpts \
    ${bwaRG} \
    $args \
    -t $task.cpus \
    \${localIndex} \
    $reads | samtools ${samtoolsCmd} -O bam -@ ${task.cpus} -o ${prefix}_\${refName}.bam -

  echo "Bwa-mem "\$(bwa 2>&1 | grep Version | cut -d" " -f2) &> versions.txt

  """
}
