/*
 * Structural variant calling with Manta
 */

process manta {
  tag "${meta.id}"
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
  def prefix = task.ext.prefix ?: "${meta.id}"
  def args = task.ext.args ?: ''
  vcftype = "${meta.status}" == "tumor" ? "tumor" : "diploid"
  inputs = "${meta.status}" == "pair" ? "--normalBam ${bam[1]} --tumorBam ${bam[0]}" : "${meta.status}" == "tumor" ? "--tumorBam ${bam}" : "--bam  ${bam}"
  """
  ${beforeScript}
  configManta.py \
    ${inputs} \
    --reference ${fasta} \
    ${args} \
    --runDir Manta

  python Manta/runWorkflow.py -m local -j ${task.cpus}

  mv Manta/results/variants/candidateSmallIndels.vcf.gz \
    Manta_${prefix}.candidateSmallIndels.vcf.gz
  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    Manta_${prefix}.candidateSmallIndels.vcf.gz.tbi
  mv Manta/results/variants/candidateSV.vcf.gz \
    Manta_${prefix}.candidateSV.vcf.gz
  mv Manta/results/variants/candidateSV.vcf.gz.tbi \
    Manta_${prefix}.candidateSV.vcf.gz.tbi
  mv Manta/results/variants/${vcftype}SV.vcf.gz \
    Manta_${prefix}.${vcftype}SV.vcf.gz
  mv Manta/results/variants/${vcftype}SV.vcf.gz.tbi \
    Manta_${prefix}.${vcftype}SV.vcf.gz.tbi
  configManta.py --version &> versions.txt 2>&1 || true
  """
}
