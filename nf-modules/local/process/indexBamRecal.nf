/*
 * Index BQSR bam
 */

 process indexBamRecal {
   tag "${meta.id}"
   label 'samtools'
   label 'minCpu'
   label 'minMem'

   input:
   tuple val(meta), path(bqsrBam)

   output:
   tuple val(meta), path(bqsrBam), path("*bam.bai"), emit:bqsrBam
   path("versions.txt"), emit: versions

   script:
   """
   samtools index ${bqsrBam}
   samtools --version &> "versions.txt" 2>&1 || true
   """
 }
