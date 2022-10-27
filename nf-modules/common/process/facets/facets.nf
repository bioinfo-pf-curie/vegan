/*
 * Compute Cellular Fraction and Copy Numbers from Tumor Sequencing
 */

process facets{
  label 'facets'
  label 'minCpu'
  label 'medMem'

  tag "${meta.tumor_id}_vs_${meta.normal_id}"

  when:
  task.ext.when == null || task.ext.when

  input:
  tuple val(meta), path(snppileupCounts)

  output:
  path("*.{txt,pdf}"), emit: facetsResultsCh
  path("versions.txt"), emit: versions

  script:
  """
  genome=\$(basename ${params.genome} "_base")

  facets.r -i ${snppileupCounts} \\
	   --name ${meta.tumor_id}_vs_${meta.normal_id} \\
	   --assembly \$genome \\
	   --normalDepth 25 --maxDepth 1000 --ampCopy 5 --hetThres 0.25
  R -e "packageVersion('facets')" > versions.txt
  """
 }
