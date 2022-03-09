/* 
 * Qualimap - quality controls on RNA-seq BAM file
 * External parameters :
 * @ params.singleEnd :	is data	single-end sequencing ?
 */

process qualimap {
  tag "${meta.id}"
  label 'qualimap'
  label 'minCpu'
  label 'medMem'

  input:
  tuple val(meta.id), path(bam), path(bai), val(stranded)
  path gtf

  output:
  path ("${meta.id}"), emit: results
  path ("versions.txt"), emit: versions

  script:
  peOpts = params.singleEnd ? '' : '-pe'
  memory = task.memory.toGiga() + "G"
  strandnessOpts = 'non-strand-specific'
  if (stranded == 'forward') {
    strandnessOpts = 'strand-specific-forward'
  } else if (stranded == 'reverse') {
    strandnessOpts = 'strand-specific-reverse'
  }
  """
  mkdir tmp
  export _JAVA_OPTIONS=-Djava.io.tmpdir=./tmp
  qualimap \\
    --java-mem-size=$memory \\
    rnaseq \\
    -bam ${bam} \\
    -gtf ${gtf} \\
    -p ${strandnessOpts} \\
    ${peOpts} \\
    -outdir ${meta.id}
  echo \$(qualimap -h 2>&1 |  grep "QualiMap" | sed -e 's/v.//') > versions.txt
  """
}
