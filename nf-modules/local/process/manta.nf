// STEP MANTA.2 - SOMATIC PAIR

process manta {
  label 'manta'
  label 'highCpu'
  label 'highMem'

  tag "${meta.tumor_id}_vs_${meta.normal_id}"

  input:
  tuple val(meta), path(bamTumor), path(baiTumor),path(bamNormal), path(baiNormal)
  path(targetBed)
  path(fasta)
  path(fastaFai)

  output:
  tuple val("Manta"), val("${meta.id}"),    val("${meta.tumor_id}_vs_${meta.normal_id}"), path("*.vcf.gz"), path("*.vcf.gz.tbi") ,emit: vcfMantaCh
  path('versions.txt') ,emit: mantaVersionCh

  when:
  task.ext.when == null || task.ext.when

  script:
  def beforeScript = task.ext.beforeScript ?: ''
  def args = task.ext.args ?: ''
  """
    ${beforeScript}
    configManta.py \
        --normalBam ${bamNormal} \
        --tumorBam ${bamTumor} \
        --reference ${fasta} \
        ${args} \
        --runDir Manta

    python Manta/runWorkflow.py -m local -j ${task.cpus}

    mv Manta/results/variants/candidateSmallIndels.vcf.gz \
        Manta_${meta.tumor_id}_vs_${meta.normal_id}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        Manta_${meta.tumor_id}_vs_${meta.normal_id}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        Manta_${meta.tumor_id}_vs_${meta.normal_id}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        Manta_${meta.tumor_id}_vs_${meta.normal_id}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/diploidSV.vcf.gz \
        Manta_${meta.tumor_id}_vs_${meta.normal_id}.diploidSV.vcf.gz
    mv Manta/results/variants/diploidSV.vcf.gz.tbi \
        Manta_${meta.tumor_id}_vs_${meta.normal_id}.diploidSV.vcf.gz.tbi
    mv Manta/results/variants/somaticSV.vcf.gz \
        Manta_${meta.tumor_id}_vs_${meta.normal_id}.somaticSV.vcf.gz
    mv Manta/results/variants/somaticSV.vcf.gz.tbi \
        Manta_${meta.tumor_id}_vs_${meta.normal_id}.somaticSV.vcf.gz.tbi
    configManta.py --version &> versions.txt 2>&1 || true
    """
}
