"""
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk)
https://github.com/cov-ert/fastafunk

This module collapses fasta sequences into consensus sequences based on groupings by some index.

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
import re
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo

from fastafunk.utils import *

def create_consensus(in_fasta,in_metadata,index_field,index_column,clade_file,out_fasta,log_file):
    metadata_dic = {}
    phylotype_dic = {}
    seq_dic = {}
    consensus_dic = {}
    output_folder = os.path.dirname(out_fasta)
    log_handle = get_log_handle(log_file, out_fasta)
    log_handle.write("Output folder: %s" %output_folder)

    with open(in_metadata,"r",encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [name.lower() for name in reader.fieldnames]
        metadata = [r for r in reader]

    for items in metadata:
        if index_field.lower() not in reader.fieldnames or index_column.lower() not in reader.fieldnames:
            sys.exit("Column name not in metadata file, please re-check metadata file and reinsert a column name.")
        else:
            metadata_dic[items[index_column]] = items[index_field.lower()]

    for record in SeqIO.parse(in_fasta, 'fasta'):
        seq_dic[record.id]= record.seq

    if os.path.isfile(clade_file):
        with open(clade_file,"r") as f:
            for line in f:
                phylotype_dic[line.rstrip()] = []
    else:
        for trait in metadata_dic.values():
            if trait not in phylotype_dic.keys():
                phylotype_dic[trait] = []

    trait_order = list(phylotype_dic.keys())
    trait_order.sort(key=lambda x: re.sub("[^A-Z0-9]", "",x),reverse=True)

    if len(set(metadata_dic.keys())&set(seq_dic.keys())) == 0:
        sys.exit("No matching sequence name with metadata name. Program Exit")

    for cluster in trait_order:
        for seq_id,phylotype in metadata_dic.items():
            cluster_type = cluster.split(".")
            cluster_length = len(cluster_type)
            phylo_type = phylotype.split(".")
            if len(phylo_type) < cluster_length:
                continue
            if phylo_type[:cluster_length] == cluster_type:
                if seq_id in seq_dic.keys():
                    phylotype_dic[cluster].append([seq_id,seq_dic[seq_id],phylotype])
                    del seq_dic[seq_id]

    for key,value in phylotype_dic.items():
        log_handle.write("%s\t%i" %(key,len(value)))

    for seq in seq_dic.keys():
        log_handle.write("Sequence " + seq + " did not find any matches to metadata file.\n")

    for key in phylotype_dic.keys():
        if len(phylotype_dic[key]) > 2:
            outfile_name = output_folder + trait + ".fasta"
            outfile = open(outfile_name,"w")
            for sequences in phylotype_dic[key]:
                record = SeqRecord(sequences[1],id=sequences[0],description="")
                SeqIO.write(record, outfile, "fasta")
            outfile.close()
            alignment_name = outfile_name[:-6] + "_alignment.fasta"
            align_command = "mafft " + outfile_name + " > " + alignment_name
            os.system(align_command)
            os.remove(outfile_name)
            alignment = AlignIO.read(alignment_name, 'fasta')
            consensus_name = key + "_consensus"
            summary_align = AlignInfo.SummaryInfo(alignment)
            consensus_seq = summary_align.dumb_consensus(threshold=0.0,ambiguous='N')
            consensus_dic[consensus_name] = consensus_seq
            os.remove(alignment_name)
        else:
            log_handle.write("Phylotype " + key + " does not have 2 or more sequences for an alignment to work.\n")
    close_handle(log_handle)

    out_handle = get_out_handle(out_fasta)
    for key, value in consensus_dic.items():
        record = SeqRecord(value,id=key,description="")
        SeqIO.write(record, out_handle, "fasta")
    close_handle(out_handle)
