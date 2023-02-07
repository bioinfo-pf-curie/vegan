/*
 * Create intervals files for parallel computing
 * From https://github.com/nf-core/sarek/blob/master/modules/local/create_intervals_bed/main.nf
 */

process createIntervalBeds {
  label 'unix'
  label 'minCpu'
  label	'minMem'

  input:
  path(intervals)

  output:
  path('*.bed'), emit: bed
  path("versions.txt"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  // If the interval file is BED format, the fifth column is interpreted to
  // contain runtime estimates, which is then used to combine short-running jobs
  if (intervals.toString().toLowerCase().endsWith("bed")) {
    """
    awk -v FS="\t" '{
      t = \$5  # runtime estimate
      if (t == "") {
        # no runtime estimate in this row, assume default value
        t = (\$3 - \$2) / ${params.nucleotidesPerSecond}
      }
      if (name == "" || (chunk > 600 && (chunk + t) > longest * 1.05)) {
        # start a new chunk
        name = sprintf("%s_%d-%d.bed", \$1, \$2+1, \$3)
        chunk = 0
        longest = 0
      }
      if (t > longest)
        longest = t
      chunk += t
      print \$0 > name
    }' ${intervals}
    echo "awk "\$(awk --version | head -1 | awk '{print \$NF}') > versions.txt
    """
  }else if (intervals.toString().toLowerCase().endsWith("interval_list")) {
    """
    grep -v '^@' ${intervals} | awk -v FS="\t" '{
      name = sprintf("%s_%d-%d", \$1, \$2, \$3);
      printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
    }'
    echo "awk "\$(awk --version | head -1 | awk '{print \$NF}') > versions.txt
    """
  }else{
    """
    awk -v FS="[:-]" '{
      name = sprintf("%s_%d-%d", \$1, \$2, \$3);
      printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
    }' ${intervals}
    echo "awk "\$(awk --version | head -1 | awk '{print \$NF}') > versions.txt
    """
  }
}
