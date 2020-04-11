"""
This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

from fastafunk.annotate import *

def run(options):

    annotate(
        options.in_fasta,
        options.in_metadata,
        options.index_column,
        options.index_field,
        options.out_fasta,
        options.out_metadata,
        options.header_delimiter,
        options.add_cov_id,
        options.log_file,
    )