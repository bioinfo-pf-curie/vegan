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
  path splan
  path metadata
  path multiqcConfig
  path ('fastqc/*')
  path ('mapping/*')
  path ('mapping/*')
  path ('metrics/*')
  path ('metrics/*')
  path ('preseq/*')
  path ('preprocessing/*')
  path ('preprocessing/*')
  path ('preprocessing/*')
  path ('coverage/*')
  path ('identito/*')
  path ('vcfMetrics/*')
  path ('tstv/*')
  path ('snpEff/*')
  path ('msi/*')
  path ('tmb/*')
  path ('softwareVersions/*')
  path ('workflowSummary/*')
  path warnings

  output:
  path splan, emit: splan
  path "*report.html", emit: report
  path "*_data", emit: data

  script:
  def prefix = task.ext.prefix ?: ""
  def args = task.ext.args ?: ""

  def metadataOpts = params.metadata ? "--metadata ${metadata}" : ""
  def designOpts= params.design ? "-d ${params.design}" : ""
  def splanOpts = params.samplePlan ? "--splan ${params.samplePlan}" : ""
  def isPE = params.singleEnd ? 0 : 1
  def modulesList = "-m custom_content -m cutadapt -m fastqc -m preseq -m picard -m gatk -m bcftools -m snpeff -m picard -m mosdepth -m samtools"
  def warn = warnings.name == 'warnings.txt' ? "--warn warnings.txt" : ""
  """
  stats2MultiQC.sh -s ${splan} ${designOpts} -p ${isPE}
  medianReadNb="\$(sort -t, -k3,3n mqc.stats | awk -F, '{a[i++]=\$3;} END{x=int((i+1)/2); if (x<(i+1)/2) printf "%.0f", (a[x-1]+a[x])/2; else printf "%.0f",a[x-1];}')"
  mqc_header.py --name "VEGAN" --version ${workflow.manifest.version} ${metadataOpts} ${splanOpts} ${warn} --nbreads \${medianReadNb} > multiqc-config-header.yaml
  multiqc . ${args} -c $multiqcConfig -c multiqc-config-header.yaml $modulesList --filename ${prefix}multiqc_report
  """    
}
