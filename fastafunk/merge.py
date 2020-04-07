from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import csv
import argparse

def merge_fasta(in_fasta, in_metadata, out_fasta):
    metadata_dictionary = {}
    sequence_dictionary = {}
    merged_file = open(out_fasta,"w")
    log_file = open(out_fasta+".log","w")
    with open(in_metadata) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            metadata_dictionary[row[0]] = row[1:]

    for fasta_file in in_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in metadata_dictionary.keys() and record.id not in sequence_dictionary.keys():
                sequence_dictionary[record.id] = record.seq
            elif record.id not in metadata_dictionary.keys():
                log_file.write(record.id + " sequence is not in metadata file or the name is wrong (in file " + fasta_file + ")\n")
            elif record.id in sequence_dictionary.keys():
                log_file.write(record.id + " is a duplicate (in file " + fasta_file + ")\n")

    for key, value in sequence_dictionary.items():
        records = SeqRecord(value, key, description= '')
        SeqIO.write(records, merged_file, "fasta")

    merged_file.close()
    log_file.close()