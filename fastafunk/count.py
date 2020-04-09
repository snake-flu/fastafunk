"""
Name: count.py
Author: Rachel Colquhoun
Date: 08 April 2020
Description: Counts the size each group, as specified by index column(s) in a metadata file.

This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import csv
import sys
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from fastafunk.utils import *

def count_groups(in_metadata,index_column,log_file):
    log_handle = get_log_handle(log_file, None)

    metadata = load_metadata(in_metadata, None, None)
    get_groups(metadata, index_column, log_handle)

    close_handle(log_handle)




