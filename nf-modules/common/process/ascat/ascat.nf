/*
* STEP ASCAT.3 - ASCAT : infers tumour purity and ploidy and calculates allele-specific copy number propaths.

*/

process ascat {
  label 'ascat'
  label 'lowMem'
  label 'minCpu'

  tag "${meta.tumor_id}_vs_${meta.normal_id}"

  input:
  tuple val(meta),
        path(bafNormal), path(logrNormal),
        path(bafTumor), path(logrTumor)
  path(acLociGC)

  output:
  tuple val("ASCAT"), val(meta), path("${meta.tumor_id}.*.{png,txt}"),emit: ascatOutCh
  path("versions.txt"),emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  purityPloidy = (params.ascatPurity && params.ascatPloidy) ? "--purity ${params.ascatPurity} --ploidy ${params.ascatPloidy}" : ""
  """
  for f in *BAF *LogR; do sed 's/chr//g' \$f > tmppath; mv tmppath \$f;done
  runAscat.r \
    --tumorbaf ${bafTumor} \
    --tumorlogr ${logrTumor} \
    --normalbaf ${bafNormal} \
    --normallogr ${logrNormal} \
    --tumorname ${meta.tumor_id} \
    --basedir ${workflow.projectDir} \
    --gcfile ${acLociGC} \
    --gender ${meta.sex} \
    ${purityPloidy}
  R -e "packageVersion('ASCAT')" > versions.txt
  """
}
