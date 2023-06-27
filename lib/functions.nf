// VEGAN : custom functions 

/**
 * Check fraction of aligned reads
 *
 * @params prefix
 * @params logs
 */

def checkAlignmentPercent(prefix, logs) {
  def percentAligned = 0;
  logs.eachLine { line ->
    if ((matcher = line =~ /([\d\.]+) \+ ([\d\.]+) mapped \s*/)) {
      nbAligned = matcher[0][1]

    } else if ((matcher = line =~ /([\d\.]+) \+ ([\d\.]+) in total \s*/)) {
      nbTotal = matcher[0][1]
    }
  }
  percentAligned = nbAligned.toFloat() / nbTotal.toFloat() * 100
  if(percentAligned.toFloat() <= '2'.toFloat() ){
      log.info "#################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($prefix)    >> ${percentAligned}% <<"
      //skippedPoorAlignment << $prefix
      return false
  } else {
      log.info "          Passed alignment > ${prefix} >> ${percentAligned}% <<"
      return true
  }
}

/**
 * Load VEGAN design file
 *
 * @params design
 */

  def loadDesign(designPath) {
    def designFile = designPath ? file(designPath) : null
    def separator = designFile.toString().endsWith(".csv") ? ',' : designFile.toString().endsWith(".tsv") ? '\t':  ''
    if (designFile) {
      return Channel.of(designFile)
        .splitCsv(sep: separator, header:true)
        .map { row -> [row.TUMOR_ID.replaceAll("[-. ]","_"), 
                       row.GERMLINE_ID.replaceAll("[-. ]","_"), 
                       row.PAIR_ID.replaceAll("[-. ]","_"), row.SEX] }
    } else {
      return Channel.empty()
    }
  }

  /**
   * Check if PON is defined for tumor only samples
   *
   */

  def checkTumorOnly(ids, params){
    if (ids.size() > 0 && !params.pon){
      exit 1, "Missing --pon option for tumor only samples"
    }
  }
