#!/usr/bin/env python
# coding: utf-8
#
#  This file is part of multiQC software. Using a metadata file, it writes part of the multiQC config file
#
#  Copyright (c) 2018 - Institut Curie
#
#  File author(s):
#      Dimitri Desvillechabrol <dimitri.desvillechabrol@curie.fr>,
#      Nicolas Servant <nicolas.servant@curie.fr>
#
#  Distributed under the terms of the CeCILL-B license.
#  The full license is in the LICENSE file, distributed with this software.
#
##############################################################################

import os
import re
import argparse
from collections import OrderedDict

ARGS = OrderedDict({
    ("-n", "--name"): {"help": "Pipeline name"},
    ("-v", "--version"): {"help": "Pipeline version", "default": ""},
    ("-m", "--metadata"): {"help": "Metatdata file", "default": None},
    ("-s", "--splan"): {"help": "Sample plan", "default": None},
    ("-x", "--nbreads"): {"help": "Number of reads to display on the graph", "default": 0}
})


def make_multiqc_header(cli_args):
    """
    Format multiqc header according to CLI args
    :param cli_args:
    :return:
    """
    splan_desc = """
sample_names_rename_buttons:
    - 'Sample ID'
    - 'Sample Name'
sample_names_rename:\
""" if cli_args.splan else ""
    dev_desc = """
report_comment: >
  This software is currently under active development and the results have been generated with a non stable version. 
  The reliability, reproducibility and the quality of the results are therefore not guaranteed.\
""" if re.match(r".*dev$", args.version) else ""
    with open(cli_args.splan, 'r') as fp:
        splan_content = '\n' + "\n".join(
            ['    - ["{}","{}"]'.format(row[0], row[1].strip(','))
             for row in [line.split('\t' if cli_args.splan.endswith("tsv") else ',') for line in fp]
             ])
    nbreads_cont = f"""
custom_plot_config:
   preseq_plot:
      xPlotLines:
         - color: '#a9a9a9'
           value: {int(args.nbreads) / 1000000:.2f}   
           dashStyle: 'LongDash'
           width: 1
           label:
              style: {{color: '#a9a9a9'}}
              text: 'Median Reads Number'
              verticalAlign: 'top'
              y: 0\
""" if int(cli_args.nbreads) > 0 else ""
    # create rims dict
    rims_dict = OrderedDict([
        ('RIMS_ID', "RIMS code"),
        ('project_name', "Project name"),
        ('project_id', "Project ID"),
        ('runs', 'Runs'),
        ('sequencer', "Sequencing setup"),
        ('biological_application', "Application type"),
        ('nature_of_material', 'Material'),
        ('protocol', 'Protocol'),
        ('bed', 'BED of targets'),
        ('technical_contact', "Main contact"),
        ('team_leader|unit', "Team leader"),
        ('ngs_contact', "Contact E-mail")
    ])

    # get data from metadata
    with open(cli_args.metadata, 'r') as fp:
        meta_dict = {row[0]: row[1].strip() for row in [line.split('\t') for line in fp]}
        # add ngs mail if no agent was set
    meta_dict['ngs_contact'] = 'ngs.lab@curie.fr'
    meta_cont = "\n".join([
        '    - {}: "{}"'.format(value, meta_dict[key])
        for key, value in rims_dict.items() if key in meta_dict])

    return f"""\
title: '{cli_args.name}'
subtitle: Institut Curie NGS/Bioinformatics core facilities
intro_text: >
  This report has been generated by the {cli_args.name} analysis pipeline ({cli_args.version})\
{dev_desc}
custom_logo: '{os.sep.join([os.path.dirname(os.path.realpath(__file__)), '../assets/institutCurieLogo.png'])}'
custom_logo_title: Institut Curie
custom_logo_url: https://science.curie.fr/plateformes/sequencage-adn-haut-debit-ngs/\
{splan_desc}{splan_content}{nbreads_cont}
report_header_info:
{meta_cont}
"""


def get_args(arg_defs):
    """
    Get args from CLI
    :param arg_defs:
    :return:
    """
    arg_parser = argparse.ArgumentParser()
    for arg_def, kwargs in arg_defs.items():
        arg_parser.add_argument(*arg_def, **kwargs)
    return arg_parser, arg_parser.parse_args()


if __name__ == '__main__':
    parser, args = get_args(ARGS)
    print(make_multiqc_header(args))
