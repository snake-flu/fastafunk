"""
Name: unwrap.py
Author: Xiaoyu Yu
Date: 10 April 2020
Description: Unwrap fasta sequences in two line format instead of line wrapping format. Program will also look for spaces and remove them within sequences which are not suppose to be there and also flag them in log file.

This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna
from fastafunk.utils import *

def unwrap_fasta(in_fasta, out_fasta, log_file):
    """
    Remove whitespace from the sequences (Undo line wrapping).

    :param in_fasta: Fasta file with sequences that needs to be filtered for white spacea or in line wrapping state. (Required)
    :param out_fasta: Output fasta file for filtered sequences with removed wrapping format
    (Default: unwrap.fasta). (Optional)
    :param log_file: Output log file (Default: stdout). (Optional)

    :return:
    """
    if not in_fasta:
        in_fasta = [""]

    out_handle = get_out_handle(out_fasta)
    log_handle = get_log_handle(log_file, out_fasta)

    for fasta_file in in_fasta:
        fasta_handle = get_in_handle(fasta_file)
        for record in SeqIO.parse(fasta_handle, "fasta"):
            print("Sequence " + record.id + " changed to fasta 2 line format",file=log_handle)
            SeqIO.write(record, out_handle, "fasta-2line")
        close_handle(fasta_handle)

    close_handle(out_handle)
    close_handle(log_handle)
