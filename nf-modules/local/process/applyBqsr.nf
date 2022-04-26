/*
 * Apply BQSR - detects systematic errors made by the sequencing machine
 */

 process applyBQSR {
   tag "${meta.id}"
   label 'gatk'
   label 'lowMem'
   label 'lowCpu'

   input:
   tuple val(meta), path(bamFiltered), path(bamFilteredBai), path(recalTable)
   path(bed)
   path(fasta)
   path(fastaFai)
   path(dict)

   output:
   tuple val(meta),path("*.recal.bam"), emit:bqsrBam
   path("versions.txt"), emit: versions

   script:
   def args = task.ext.args ?: ''
   def prefix = task.ext.prefix ?: "${meta.id}"
   """
   gatk --java-options -Xmx${task.memory.toGiga()}g \
       ApplyBQSR \
       -R ${fasta} \
       --input ${bamFiltered} \
       --output ${prefix}.recal.bam \
       ${args} \
       --bqsr-recal-file ${recalTable}
   gatk ApplyBQSR --help &> "versions.txt" 2>&1 || true
   """
 }
