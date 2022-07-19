/*
 * Structural variant calling with Manta
 */

process manta {
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
  vcftype = "${meta.status}" == "tumor" ? "tumor" : "diploid"
  fileID = "${meta.status}" == "pair" ? "${meta.tumor_id}_vs_${meta.normal_id}" : "${meta.status}" == "tumor" ? "${meta.tumor_id}" : "${meta.normal_id}"
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
    Manta_${fileID}.candidateSmallIndels.vcf.gz
  mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
    Manta_${fileID}.candidateSmallIndels.vcf.gz.tbi
  mv Manta/results/variants/candidateSV.vcf.gz \
    Manta_${fileID}.candidateSV.vcf.gz
  mv Manta/results/variants/candidateSV.vcf.gz.tbi \
    Manta_${fileID}.candidateSV.vcf.gz.tbi
  mv Manta/results/variants/${vcftype}SV.vcf.gz \
    Manta_${fileID}.${vcftype}SV.vcf.gz
  mv Manta/results/variants/${vcftype}SV.vcf.gz.tbi \
    Manta_${fileID}.${vcftype}SV.vcf.gz.tbi
  configManta.py --version &> versions.txt 2>&1 || true
  """
}
