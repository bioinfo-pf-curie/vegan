/*
* STEP ASCAT.1 - ALLELECOUNTER
*/

process alleleCounter {
  label 'ascat'
  label 'lowMem'
  label 'minCpu'

  tag "${fileID}"

  input:
  tuple val(meta), path(bam), path(bai)
  path(acLoci)
  path(dict)
  path(fasta)
  path(fastaFai)

  output:
  tuple val(meta), path("${fileID}.alleleCount"), emit: alleleCounterOutCh
  path("versions.txt"),emit: version

  when:
  task.ext.when == null || task.ext.when

  script:
  fileID = "${meta.status}" == "tumor" ? "${meta.tumor_id}" : "${meta.normal_id}"
  """
  alleleCounter \
    -l ${acLoci} \
    -r ${fasta} \
    -b ${bam} \
    -o ${fileID}.alleleCount;
  alleleCounter --version &> versions.txt 2>&1 || true
  """
}
