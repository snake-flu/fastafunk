
from fastafunk.add_columns import *

def run(options):

    add_columns(options.in_metadata,
                options.in_data,
                options.index_column,
                options.join_on,
                options.new_columns,
                options.out_metadata,
                options.where_column,
                options.log_file
                            )
