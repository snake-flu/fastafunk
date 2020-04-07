"""
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk)
https://github.com/cov-ert/fastafunk

This module extracts sequences based on matches to the metadata.

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

def extract_fasta(in_fasta, in_metadata, out_fasta):
    metadata_dictionary = {}

    with open(in_metadata) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            metadata_dictionary[row[0]] = row[1:]

    extract_file = open(out_fasta,"w")
    log_file = open(out_fasta+".log","w")

    for record in SeqIO.parse(in_fasta, "fasta"):
        if record.id in metadata_dictionary.keys():
            SeqIO.write(record, extract_file, "fasta")
        else:
            log_file.write("Sequence " + record.id + " removed due to no match to metadata\n")
    extract_file.close()
    log_file.close()