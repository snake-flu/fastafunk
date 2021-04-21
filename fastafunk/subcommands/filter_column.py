
from fastafunk.filter_column import *

def run(options):

    filter_column(options.in_metadata,
                options.column,
                options.out_metadata,
                options.is_true,
                options.is_false,
                options.log_file)
