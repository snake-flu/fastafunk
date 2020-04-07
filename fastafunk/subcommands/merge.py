from fastafunk.merge import *

def run(options):

    merge_fasta(
        options.in_fasta,
        options.in_metadata,
        options.out_fasta,
        options.log_file
    )