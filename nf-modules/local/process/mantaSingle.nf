/*
 * Structural variant calling with Manta
 */

process mantaSingle {
  tag "${meta.status}"
  label 'manta'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(meta), path(bam), path(bai)
  path(targetBed)
  path(fasta)
  path(fastaFai)

  output:
  tuple val(meta), val("Manta"), path("*.vcf.gz"), path("*.vcf.gz.tbi"), emit:svVcf
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def beforeScript = task.ext.beforeScript ?: ''
  def args = task.ext.args ?: ''

  inputbam = "${meta.status}" == "tumor" ? "--tumorBam" : "--bam"
  vcftype = "${meta.status}" == "tumor" ? "tumor" : "diploid"
  """
  echo "${meta.status}"

  ${beforeScript}
  configManta.py \
    ${inputbam} ${bam} \
    --reference ${fasta} \
    ${args} \
    --runDir Manta

  python Manta/runWorkflow.py -m local -j ${task.cpus}

  mv Manta/results/variants/candidateSmallIndels.vcf.gz \
    Manta_${meta.id}.candidateSmallIndels.vcf.gz
  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    Manta_${meta.id}.candidateSmallIndels.vcf.gz.tbi
  mv Manta/results/variants/candidateSV.vcf.gz \
    Manta_${meta.id}.candidateSV.vcf.gz
  mv Manta/results/variants/candidateSV.vcf.gz.tbi \
    Manta_${meta.id}.candidateSV.vcf.gz.tbi
  mv Manta/results/variants/${vcftype}SV.vcf.gz \
    Manta_${meta.id}.${vcftype}SV.vcf.gz
  mv Manta/results/variants/${vcftype}SV.vcf.gz.tbi \
    Manta_${meta.id}.${vcftype}SV.vcf.gz.tbi
  configManta.py --version &> versions.txt 2>&1 || true
  """
}
