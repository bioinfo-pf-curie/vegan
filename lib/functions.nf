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
    //def designExt = designPath ? getExtension(designPath, ["tsv", "csv"]) : ""
    //def separator = (designExt == 'tsv') ? '\t' : (designExt == 'csv') ? ',' : ''
    def separator = designFile.toString().endsWith(".csv") ? ',' : designFile.toString().endsWith(".tsv") ? '\t':  ''
    if (designFile) {
      return Channel.of(designFile)
        .splitCsv(sep: separator, header: ['germlineId', 'tumorId', 'pairId', 'sex'])
        .map { row -> [row.germlineId, row.tumorId, row.pairId, row.sex] }
    } else {
      return Channel.empty()
    }
  }
