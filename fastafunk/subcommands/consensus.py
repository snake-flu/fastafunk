from fastafunk.consensus import *

def run(options):

    create_consensus(options.in_fasta,
                     options.in_metadata,
                     options.index_field,
                     options.index_column,
                     options.clade_file,
                     options.out_fasta
                     )