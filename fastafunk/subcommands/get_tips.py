"""
This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

from fastafunk.get_tips import *

def run(options):

    get_tips(
        options.in_metadata,
        options.in_tree,
        options.out_tips,
        options.low_memory
    )
