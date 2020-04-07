from fastafunk.remove import *

def run(options):

    remove_fasta(
        options.in_fasta,
        options.in_metadata,
        options.out_fasta,
        options.log_file
    )