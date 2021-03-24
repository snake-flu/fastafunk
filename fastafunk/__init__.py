"""
This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

from pkg_resources import get_distribution

try:
    __version__ = get_distribution("fastafunk").version
except:
    __version__ = "local"

__all__ = ["consensus", "extract", "merge", "remove", "split", "count", "subsample", "annotate",
           "unwrap","strip", "new", "add_columns", "fetch", "shuffle", "drop_columns"]

from fastafunk import *
