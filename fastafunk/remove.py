"""
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk)
https://github.com/cov-ert/fastafunk

This module removes sequences based on matches to the metadata.

This file is part of Fastafunk. Fastafunk is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Fastafunk is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Fastafunk. If
not, see <http://www.gnu.org/licenses/>.
"""

import csv
from Bio import SeqIO

from fastafunk.utils import *

def remove_fasta(in_fasta, in_metadata, out_fasta, log_file):
    if not in_fasta:
        in_fasta = [""]
    metadata_dictionary = metadata_to_dict(in_metadata)

    out_handle = get_out_handle(out_fasta)
    log_handle = get_log_handle(log_file, out_fasta)

    for fasta_file in in_fasta:
        fasta_handle = get_in_handle(fasta_file)
        for record in SeqIO.parse(fasta_handle, "fasta"):
            if record.id not in metadata_dictionary.keys():
                SeqIO.write(record, out_handle, "fasta")
            else:
                print("Sequence " + record.id + " removed due to match to metadata", file=log_handle)
        close_handle(fasta_handle)

    close_handle(out_handle)
    close_handle(log_handle)