#!/usr/bin/env python
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="version file", type=str, default='')
args = parser.parse_args()

versions = {}
with open(args.input) as f:
    for line in f:
        if line.strip():
            (key, val) = line.strip().split()
            if key in versions.keys():
                val_list = versions[str(key)]
                if val not in val_list:
                    val_list.append(val)
                    versions[str(key)] = val_list
            else:
                versions[str(key)] = [val]
# Dump to YAML
print ('''
id: 'software_versions'
section_name: 'Software Versions'
section_href: 'https://gitlab.curie.fr/rnaseq'
plot_type: 'html'
description: 'are collected at run time from the software output.'
data: |
    <dl class="dl-horizontal">
''')
for k,v in versions.items():
    print("        <dt>{}</dt><dd><samp>{}</samp></dd>".format(k,','.join(v)))
print ("    </dl>")
