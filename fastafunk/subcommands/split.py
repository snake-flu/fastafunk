from fastafunk.split import *

def run(options):

    split_fasta(
        options.in_fasta,
        options.in_metadata,
        options.index_field,
        options.index_column,
        options.out_folder
    )
