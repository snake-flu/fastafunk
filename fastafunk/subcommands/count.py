from fastafunk.count import *

def run(options):

    count_groups(options.in_metadata,
          options.index_column,
          options.log_file
         )