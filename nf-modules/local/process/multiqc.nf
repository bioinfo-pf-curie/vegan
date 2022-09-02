/*
 * MultiQC for RNA-seq report
 * External parameters :
 * @ params.singleEnd :	is data	single-end sequencing ?
 */

process multiqc {
  label 'multiqc'
  label 'minCpu'
  label 'lowMem'

  input:
  val customRunName
  path splan
  path metadata
  path multiqcConfig
  path ('fastqc/*')
  path ('mapping/*')
  path ('mapping/*')
  path ('preseq/*')
  path ('preprocessing/*')
  path ('preprocessing/*')
  path ('coverage/*')
  path ('coverage/*')
  path ('identito/*')
  path ('vcfMetrics/*')
  path ('vcfMetrics/*')
  path ('softwareVersions/*')
  path ('workflowSummary/*')
  path warnings

  output:
  path splan, emit: splan
  path "*report.html", emit: report
  path "*_data", emit: data

  script:
  rtitle = customRunName ? "--title \"$customRunName\"" : ''
  rfilename = customRunName ? "--filename " + customRunName + "_rnaseq_report" : "--filename rnaseq_report"
  metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
  designOpts= params.design ? "-d ${params.design}" : ""
  splanOpts = params.samplePlan ? "--splan ${params.samplePlan}" : ""
  isPE = params.singleEnd ? 0 : 1
    
  modulesList = "-m custom_content -m fastqc -m preseq -m picard -m gatk -m bcftools -m snpeff -m picard -m mosdepth"
  warn = warnings.name == 'warnings.txt' ? "--warn warnings.txt" : ""
  """
  apStats2MultiQC.sh -s ${splan} ${designOpts} ${isPE}
  medianReadNb="\$(sort -t, -k3,3n mqc.stats | awk -F, '{a[i++]=\$3;} END{x=int((i+1)/2); if (x<(i+1)/2) printf "%.0f", (a[x-1]+a[x])/2; else printf "%.0f",a[x-1];}')"
  mqc_header.py --name "VEGAN" --version ${workflow.manifest.version} ${metadataOpts} ${splanOpts} ${warn} --nbreads \${medianReadNb} > multiqc-config-header.yaml
  multiqc . -f $rtitle $rfilename -c $multiqcConfig -c multiqc-config-header.yaml $modulesList
  """    
}
