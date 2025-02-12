/*
 * Facets flow
 */

include { facetsPileup } from '../../common/process/facets/facetsPileup'
include { facets} from '../../common/process/facets/facets'
include { annotateCNV } from '../../local/process/annotateCNV'

workflow facetsFlow {

    take:
    bam
    dbsnp
    gtf

    main:
    chVersions = Channel.empty()

    facetsPileup(
      bam,
      dbsnp
    )

    facets(
      facetsPileup.out.snppileup
    )
    chVersions = chVersions.mix(facets.out.versions)

    annotateCNV(
      facets.out.segmentsTransformed,
      gtf
    )

    emit:
    segmentsTransformed = facets.out.segmentsTransformed
    segmentsAnnotated = annotateCNV.out.cnv
    versions = chVersions
}