"""
Name: strip.py
Author: Xiaoyu Yu
Date: 10 April 2020
Description: Strip sites within sequences based on various options set by user. 

Options include:
    --gap (Remove "-")
    --ambiguity (Remove "N")
    --missing (Remove "?")
    --keep_alignment (Remove all gaps if it is shared between every sequence within the alignment)
    --front (Strip only at the front of the sequence)
    --back (Strip only at the end of the sequence)

This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import re
import numpy as np
from functools import reduce
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from fastafunk.utils import *


def strip_gap_alignment(in_fasta,keep_alignment):
    gap_position_list = []
    for record in SeqIO.parse(in_fasta, "fasta"):
        sequence = str(record.seq)
        if sequence.count("-") == 0:
            gap_position_list = []
            break
        gap_position_list.append([m.start() for m in re.finditer(sequence,"-")])    
    #TODO concatenating lsit of list for indexes. Joining characters based on missing index

#TODO fix functions
def strip_gap(sequence,seq_id,orientation,log_file):
    strip_seq = sequence.strip("-")
    print("Sequence " + seq_id + " contains gaps and has been removed",file=log_file)
    return strip_seq

def strip_ambiguity(sequence,seq_id,orientation,log_file):
    strip_seq = sequence.strip("N")
    print("Sequence " + seq_id + " contains ambiguous sites and has been removed",file=log_file)
    return strip_seq

def strip_missing(sequence,seq_id,orientation,log_file):
    strip_seq = sequence.strip("?")
    print("Sequence " + seq_id + " contains missing sites and has been removed",file=log_file)
    return strip_seq

def strip_fasta(in_fasta,gap,ambiguity,missing,keep_alignment,front,back,out_fasta,log_file):
    """
    Strip sites within sequences based on various options set by user. 

    :param in_fasta: Fasta file with sequences that needs to be splitted according to criteria to create a consensus
    :param gap: Remove all gaps within sequence ("-")
    :param ambiguity: Remove all ambiguious sites within sequence ("N")
    :param missing: Remove all missing sites within sequence ("?")
    :param keep_alignment: Remove only gaps ("-")
    :param front: Remove only the front portion of the sequence ("-")
    :param back: Remove only the back portion of the sequence ("-")
    :param out_fasta: Output fasta file with consensus sequences of all groups based on trait (Default: consensus.fasta). (Optional)
    :param log_file: Output log file (Default: stdout). (Optional)

    :return:
    """
    if not in_fasta:
        in_fasta = [""]

    out_handle = get_out_handle(out_fasta)
    log_handle = get_log_handle(log_file, out_fasta)

    for fasta_file in in_fasta:
        sequence_dictionary = {}
        fasta_handle = get_in_handle(fasta_file)
        for record in SeqIO.parse(fasta_handle, "fasta"):
            strip_seq = record.seq
            #TODO function and orientation list function loop
            sequence_dictionary[record.id] = strip_seq            
        close_handle(fasta_handle)

        for seq_id,seq in sequence_dictionary.items():
            record = SeqRecord(seq,id=seq_id,description="")
            SeqIO.write(record, out_handle, "fasta-2line")

    close_handle(out_handle)
    close_handle(log_handle)