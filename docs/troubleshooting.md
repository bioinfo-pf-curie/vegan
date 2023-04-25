# Troubleshooting

## Facets in tumor only mode

`Facets` requires a normal to be run. However, the normal can be unmatch.  
In this case, the easier would be to start from bam files, with a design file specifying 
unmatch pairs, and to run VEGAN with the option `--facetsOpts '--unmatch --hetThres 0.1'`  
In this mode, heterogygous SNPs are called using tumor reads only and logOR calculations are different.
That's why it is recommanded to descrease the minimal VAF value to call a SNP heterozygous (`--hetThres`)

