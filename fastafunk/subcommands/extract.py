from fastafunk.extract import *

def run(options):

    extract_fasta(
        options.in_fasta,
        options.in_metadata,
        options.out_fasta,
        options.log_file
    )