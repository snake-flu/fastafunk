import csv
from Bio import SeqIO
import argparse

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