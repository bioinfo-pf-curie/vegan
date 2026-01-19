/*
 * Structural variant calling with Manta
 */

process manta {
  tag "${meta.id}"
  label 'manta'
  label 'highCpu'
  label 'highMem'

  input:
  tuple val(meta), path(tumor_bam), path(tumor_bai), path(normal_bam), path(normal_bai)
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
  vcftype = "${meta.status}" == "tumor" ? "" : "diploid"
  """
  ${beforeScript}
  configManta.py \
    --normalBam ${normal_bam} \
    --tumorBam ${tumor_bam} \
    --reference ${fasta} \
    ${args} \
    --runDir Manta

  python Manta/runWorkflow.py -m local -j ${task.cpus}

  mv Manta/results/variants/somaticSV.vcf.gz \
  Manta_${meta.tumor_id}_vs_${meta.normal_id}.somaticSV.vcf.gz
  mv Manta/results/variants/somaticSV.vcf.gz.tbi \
  Manta_${meta.tumor_id}_vs_${meta.normal_id}.somaticSV.vcf.gz.tbi
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
