"""
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk)
https://github.com/cov-ert/fastafunk

This module merges fasta files while avoiding duplicates specified in metadata.

This file is part of Fastafunk. Fastafunk is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Fastafunk is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Fastafunk. If
not, see <http://www.gnu.org/licenses/>.
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import csv

from fastafunk.utils import *

def merge_fasta(in_fasta, in_metadata, out_fasta, log_file):
    metadata_dictionary = {}
    sequence_dictionary = {}
    out_handle = get_out_handle(out_fasta)
    log_handle = get_log_handle(log_file, out_fasta)
    with open(in_metadata) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            metadata_dictionary[row[0]] = row[1:]

    for fasta_file in in_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in metadata_dictionary.keys() and record.id not in sequence_dictionary.keys():
                sequence_dictionary[record.id] = record.seq
            elif record.id not in metadata_dictionary.keys():
                log_handle.write(record.id + " sequence is not in metadata file or the name is wrong (in file " + fasta_file + ")\n")
            elif record.id in sequence_dictionary.keys():
                log_handle.write(record.id + " is a duplicate (in file " + fasta_file + ")\n")

    for key, value in sequence_dictionary.items():
        records = SeqRecord(value, key, description= '')
        SeqIO.write(records, out_handle, "fasta")

    close_handle(out_handle)
    close_handle(log_handle)