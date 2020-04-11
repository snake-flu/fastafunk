"""
This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

from fastafunk.strip import *

def run(options):

    strip_fasta(
        options.in_fasta,
        options.gap,
        options.ambiguity,
        options.missing,
        options.keep_alignment,
        options.front,
        options.back,
        options.out_fasta,
        options.log_file
    )