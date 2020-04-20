"""
This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

from fastafunk.fetch import *

def run(options):

    fetch_fasta(
        options.in_fasta,
        options.in_metadata,
        options.index_column,
        options.filter_column,
        options.where_column,
        options.restrict,
        options.out_fasta,
        options.out_metadata,
        options.log_file,
        options.header_delimiter
    )