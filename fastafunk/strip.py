"""
Name: strip.py
Author: Xiaoyu Yu
Date: 13 April 2020
Description: Strip sites within sequences based on various options set by user. 

Options include:
    --gap (Remove "-")
    --ambiguity (Remove "N")
    --missing (Remove "?")
    --keep_alignment (Remove all gaps if it is shared between every sequence within the alignment, other functions does
    not work with this options)
    --front (Strip only at the front of the sequence)
    --back (Strip only at the end of the sequence)

Note: front and back options can be used at the same time. If both are false, function will be applied to the whole sequence

This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import re
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from fastafunk.utils import *


def strip_gap_alignment(in_fasta,out_handle,log_file):
    gap_position_dic = {}
    for record in SeqIO.parse(in_fasta, "fasta"):
        sequence = str(record.seq)
        if sequence.count("-") == 0:
            gap_position_dic = {}
            break
        string_array = np.array(list(sequence))
        gap_indexes = np.where(string_array == '-')
        gap_position_dic[record.id] = [sequence,gap_indexes[0].tolist()]

    if gap_position_dic != {}:
        gap_position_list = [value[1] for value in gap_position_dic.values()]
        common_gap = list(set.intersection(*map(set, gap_position_list)))
        common_gap.sort()
        for value in gap_position_dic.values():
            filtered_seq = "".join([char for idx, char in enumerate(value[0]) if idx not in common_gap])
            value[0] = filtered_seq

        for seq_id,seq in gap_position_dic.items():
            record = SeqRecord(Seq(seq[0]),id=seq_id,description="")
            print("Gaps at positions: " + ','.join(str(e) for e in common_gap) + " has been removed for sequence "+ seq_id,file=log_file)
            SeqIO.write(record, out_handle, "fasta-2line")
    else:
        print("There is no common gaps within the alignment to remove",file=log_file)




def strip_gap(sequence,seq_id,orientation,log_file):
    if orientation[0] == True:
        strip_seq = sequence.lstrip("-")
        print("Gaps at the front of sequence " + seq_id + " has been removed",file=log_file)
    if orientation[1] == True:
        strip_seq = sequence.rstrip("-")
        print("Gaps at the back of sequence " + seq_id + " has been removed",file=log_file)
    if orientation[0] == False and orientation[1] == False:
        strip_seq = sequence.replace("-","")
        print("Gaps in sequence " + seq_id + " has been removed",file=log_file)
    return strip_seq

def strip_ambiguity(sequence,seq_id,orientation,log_file):
    if orientation[0] == True:
        strip_seq = sequence.lstrip("N")
        print("Ambiguities (N) at the front of sequence " + seq_id + " has been removed",file=log_file)
    if orientation[1] == True:
        strip_seq = sequence.rstrip("N")
        print("Ambiguities (N) at the back of sequence " + seq_id + " has been removed",file=log_file)
    if orientation[0] == False and orientation[1] == False:
        strip_seq = sequence.replace("N","")
        print("Ambiguities (N) in sequence " + seq_id + " has been removed",file=log_file)
    return strip_seq

def strip_missing(sequence,seq_id,orientation,log_file):
    if orientation[0] == True:
        strip_seq = sequence.lstrip("?")
        print("Missing sites (?) at the front of sequence " + seq_id + " has been removed",file=log_file)
    if orientation[1] == True:
        strip_seq = sequence.rstrip("?")
        print("Missing sites (?) at the back of sequence " + seq_id + " has been removed",file=log_file)
    if orientation[0] == False and orientation[1] == False:
        strip_seq = sequence.replace("?","")
        print("Missing sites (?) in sequence " + seq_id + " has been removed",file=log_file)
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

    orientation = [front,back]
    function = [gap,ambiguity,missing]

    for fasta_file in in_fasta:
        sequence_dictionary = {}
        fasta_handle = get_in_handle(fasta_file)
        if keep_alignment == True:
            strip_gap_alignment(fasta_handle,out_handle,log_handle)
        elif not any(function):
            print("No options chosen for split function to apply on sequences", file=log_file)
        else:
            for record in SeqIO.parse(fasta_handle, "fasta"):
                sequence = str(record.seq)
                if function[0] == True:
                    sequence = strip_gap(sequence,record.id,orientation,log_handle)
                if function[1] == True:
                    sequence = strip_ambiguity(sequence,record.id,orientation,log_handle)
                if function[2] == True:
                    sequence = strip_missing(sequence,record.id,orientation,log_handle)                        
                sequence_dictionary[record.id] = sequence
            close_handle(fasta_handle)

            for seq_id,seq in sequence_dictionary.items():
                record = SeqRecord(Seq(seq),id=seq_id,description="")
                SeqIO.write(record, out_handle, "fasta-2line")

    close_handle(out_handle)
    close_handle(log_handle)