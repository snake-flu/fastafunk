import os
import csv
import re
import argparse
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo

def create_consensus(in_fasta,in_metadata,index_field,index_column,clade_file,out_fasta):
    metadata_dic = {}
    phylotype_dic = {}
    seq_dic = {}
    consensus_dic = {}
    output_folder = os.path.dirname(out_fasta)
    log_file = open(out_fasta+".log","w")
    print(output_folder)

    with open(in_metadata,"r",encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        reader.fieldnames = [name.lower() for name in reader.fieldnames]
        metadata = [r for r in reader]

    for items in metadata:
        if index_field.lower() not in reader.fieldnames or index_column.lower() not in reader.fieldnames:
            print("Column name not in metadata file, please re-check metadata file and reinsert a column name.")
            sys.exit()
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
        print("No matching sequence name with metadata name. Program Exit")
        sys.exit()

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
        print(key,len(value))

    for seq in seq_dic.keys():
        log_file.write("Sequence " + seq + " did not find any matches to metadata file.\n")

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
            log_file.write("Phylotype " + key + " does not have 2 or more sequences for an alignment to work.\n")
    log_file.close()

    consensus_file = open(out_fasta,"w")
    for key, value in consensus_dic.items():
        record = SeqRecord(value,id=key,description="")
        SeqIO.write(record, consensus_file, "fasta")
    consensus_file.close()

consensus = argparse.ArgumentParser(description='Removes sequences based on matches to the metadata')
consensus.add_argument("--in-fasta", required=True, default=None,type=str, help="<filename> [<filename> ...] one fasta files of sequences")
consensus.add_argument("--in-metadata", required=True, default=None,type=str, help="<filename> [<filename> ...] one CSV table of metadata")
consensus.add_argument("--index-field", required=True, default=None,type=str, help="<column> [<column> ...] the field(s) in the header to match the metadata")
consensus.add_argument("--index-column", required=False, default="header",type=str, help="<column> [<column> ...] the column(s) in the metadata file to use to match to the sequence")
consensus.add_argument("--clade-file", required=False, default="",type=str, help="<filename> [<filename> ...] text file including specific traits to collapse by")
consensus.add_argument("--out-fasta", required=False, default="consensus.fasta",type=str, help="<filename> output a consensus fasta file")
args = consensus.parse_args()

print(args)

in_fasta = args.in_fasta
in_metadata = args.in_metadata
out_fasta  = args.out_fasta
index_field = args.index_field
index_column = args.index_column
clade_file = args.clade_file

create_consensus(in_fasta,in_metadata,index_field,index_column,clade_file,out_fasta)