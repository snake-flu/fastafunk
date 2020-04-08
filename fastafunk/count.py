"""
Copyright 2020 Rachel Colquhoun (rachel.colquhoun@ed.ac.uk) & Xiaoyu Yu (xiaoyu.yu@ed.ac.uk)
https://github.com/cov-ert/fastafunk

This module subsamples a fasta based on counts in groupings defined by a set of index fields in metadata.

This file is part of Fastafunk. Fastafunk is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Fastafunk is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Fastafunk. If
not, see <http://www.gnu.org/licenses/>.
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




