
from fastafunk.add_columns import *

def run(options):

    drop_columns(options.in_metadata,
                options.drop_columns,
                options.out_metadata,
                options.log_file)
