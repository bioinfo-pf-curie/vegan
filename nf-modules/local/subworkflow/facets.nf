/*
 * Facets flow
 */

include { facetsPileup } from '../../common/process/facets/facetsPileup'
include { facets} from '../../common/process/facets/facets'

workflow facetsFlow {

    take:
    bam
    dbsnp

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

    emit:

    versions = chVersions
}
