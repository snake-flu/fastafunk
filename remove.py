import csv
from Bio import SeqIO
import argparse

def remove_fasta(in_fasta, in_metadata, out_fasta):
    metadata_dictionary = {}

    with open(in_metadata) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            metadata_dictionary[row[0]] = row[1:]

    filtered_file = open(out_fasta,"w")
    log_file = open(out_fasta+".log","w")

    for record in SeqIO.parse(in_fasta, "fasta"):
        if record.id not in metadata_dictionary.keys():
            SeqIO.write(record, filtered_file, "fasta-2line")
        else:
            log_file.write("Sequence " + record.id + " removed due to match to metadata\n")
    filtered_file.close()
    log_file.close()

remove = argparse.ArgumentParser(description='Removes sequences based on matches to the metadata')
remove.add_argument("--in-fasta", required=True, default=None,type=str, help="<filename> [<filename> ...] one fasta files of sequences")
remove.add_argument("--in-metadata", required=True, default=None,type=str, help="<filename> [<filename> ...] one CSV table of metadata")
remove.add_argument("--out-fasta", required=False, default="remove_by_metadata.fasta",type=str, help="<filename> output a fasta file")
args = remove.parse_args()

in_fasta = args.in_fasta
in_metadata = args.in_metadata
out_fasta  = args.out_fasta

remove_fasta(in_fasta, in_metadata, out_fasta)