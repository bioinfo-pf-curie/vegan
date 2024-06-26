/*
* STEP ASCAT.1 - ALLELECOUNTER
*/

process alleleCounter {
  label 'ascat'
  label 'lowMem'
  label 'minCpu'

  tag "${meta.id}"

  input:
  tuple val(meta), path(bam), path(bai)
  path(acLoci)
  path(dict)
  path(fasta)
  path(fastaFai)

  output:
  tuple val(meta), path("${meta.id}.alleleCount"), emit: alleleCounterOutCh
  path("versions.txt"),emit: version

  when:
  task.ext.when == null || task.ext.when

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  alleleCounter \
    -l ${acLoci} \
    -r ${fasta} \
    -b ${bam} \
    -o ${meta.id}.alleleCount;
  alleleCounter --version &> versions.txt 2>&1 || true
  echo "test"
  """
}
