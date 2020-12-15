#!/usr/bin/env python
from collections import OrderedDict
import re
import os

regexes = {
    'AlleleCount': ['v_allelecount.txt', r"(\S+)"],
    'ASCAT': ['v_ascat.txt', r"Version:       (\S+)"],
    'bcftools': ['v_bcftools.txt', r"bcftools (\S+)"],
    'BWA': ['v_bwa.txt', r"Version: (\S+)"],
    'FastQC': ['v_fastqc.txt', r"FastQC v(\S+)"],
    'FreeBayes': ['v_freebayes.txt', r"version:  v(\d\.\d\.\d+)"],
    'GATK': ['v_gatk.txt', r"Version:(\S+)"],
    'htslib': ['v_samtools.txt', r"htslib (\S+)"],
    'Manta': ['v_manta.txt', r"([0-9.]+)"],
    'MultiQC': ['v_multiqc.txt', r"multiqc, version (\S+)"],
    'Nextflow': ['v_nextflow.txt', r"(\S+)"],
    'pipeline': ['v_pipeline.txt', r"(\S+)"],
    'Qualimap': ['v_qualimap.txt', r"QualiMap v.(\S+)"],
    'R': ['v_r.txt', r"R version (\S+)"],
    'samtools': ['v_samtools.txt', r"samtools (\S+)"],
    'SnpEff': ['v_snpeff.txt', r"version SnpEff (\S+)"],
    'TIDDIT': ['v_tiddit.txt', r"TIDDIT-(\S+)"],
}
results = OrderedDict({
    key: '<span style="color:#999999;\">N/A</span>'
    for key in regexes.keys()})

# Search each file using its regex
for k, v in regexes.items():
    if os.path.exists(v[0]):
        with open(v[0]) as x:
            versions = x.read()
            match = re.search(v[1], versions)
            if match:
                results[k] = "v{}".format(match.group(1))

# Remove software set to false in results
for k in results:
    if not results[k]:
        del(results[k])

# Dump to YAML
yaml_output = '''
id: 'software_versions'
section_name: 'Software versions'
section_href: 'https://github.com/nf-core/sarek'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
'''

for k, v in results.items():
    yaml_output += "        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k, v)
print(yaml_output + "    </dl>")

# Write out regexes as csv file:
with open('software_versions.csv', 'w') as f:
    for k, v in results.items():
        f.write("{}\t{}\n".format(k, v))
