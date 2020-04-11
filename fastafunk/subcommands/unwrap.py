"""
This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

from fastafunk.unwrap import *

def run(options):

    unwrap_fasta(
        options.in_fasta,
        options.out_fasta,
        options.log_file
    )