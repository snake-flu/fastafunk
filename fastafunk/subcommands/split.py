"""
This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

from fastafunk.split import *

def run(options):

    split_fasta(
        options.in_fasta,
        options.in_metadata,
        options.index_field,
        options.index_column,
        options.lineage,
        options.lineage_csv,
        options.aliases,
        options.out_prefix,
        options.log_file
    )
