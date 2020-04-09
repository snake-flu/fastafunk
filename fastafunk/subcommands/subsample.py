"""
This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

from fastafunk.subsample import *

def run(options):

    subsample_fasta(
        options.in_fasta,
        options.in_metadata,
        options.index_field,
        options.index_column,
        options.out_fasta,
        options.out_metadata,
        options.log_file,
        options.sample_size,
        options.target_file,
        options.exclude_uk
    )