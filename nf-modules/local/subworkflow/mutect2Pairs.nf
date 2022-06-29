/*
 * Mutect2 Tumor/Normal
 */

include { mutect2 } from '../../common/process/gatk/mutect2'

workflow mutect2PairsFlow {

  take:
  bam
  bed
  fasta
  fai
  dict
  germlineResource
  germlineResourceIndex
  panelsOfNormals
  panelsOfNormalsIndex

  main:
  chVersions = Channel.empty()

  // Mutect2 tumor/bam inputs
  //[meta][tumor_bam, normal_bam],[tumor_bai, normal_bai]
  bam
    .map{ it -> [it[0], [it[1], it[3]], [it[2], it[4]]] }
    .set{ chBamMutect2 }

  chBamMutect2.view()

  mutect2(
    chBamMutect2,
    bed,
    fasta,
    fai,
    dict,
    germlineResource,
    germlineResourceIndex,
    panelsOfNormals,
    panelsOfNormalsIndex
  )

  emit:
  versions = chVersions
}
