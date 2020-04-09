"""
This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

from fastafunk.consensus import *

def run(options):

    create_consensus(options.in_fasta,
                     options.in_metadata,
                     options.index_field,
                     options.index_column,
                     options.clade_file,
                     options.out_fasta,
                     options.log_file
                     )