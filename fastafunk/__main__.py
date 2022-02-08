"""
This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import argparse
import sys

import fastafunk
import fastafunk.subcommands

def main(args=None):
    parser = argparse.ArgumentParser(
        prog="fastafunk",
        description="Miscellaneous fasta manipulation tools",
    )

    parser.add_argument("--version", action="version", version=fastafunk.__version__)
    subparsers = parser.add_subparsers(
        title="Available subcommands", help="", metavar=""
    )

    # _______________________________   common  _________________________________#
    common = argparse.ArgumentParser(add_help=False)
    # common.add_argument(
    #     '--in-metadata', dest = 'in_metadata', nargs ='+', metavar='<filename>', required = False,
    #     help='One or more CSV or TSV tables of metadata'
    # )
    # common.add_argument(
    #     '--in-fasta', dest='in_fasta', nargs='+', metavar='<filename>', required=False,
    #     help='One or more FASTA files of sequences (else reads from stdin)'
    # )
    # common.add_argument(
    #     '--index-column', dest='index_column', nargs='+', metavar='<column>', required=False,
    #     help='Column(s) in the metadata file to use to match to the sequence'
    # )
    # common.add_argument(
    #     '--index-field', dest='index_field', nargs='+', metavar='<field>', required=False,
    #     help='Field(s) in the fasta header to match the metadata (else matches column names)'
    # )
    # common.add_argument(
    #     '--where-column', dest='where_column', nargs='+', metavar='<column>=<regex>', required=False,
    #     help='Additional matches to metadata columns'
    # )
    # common.add_argument(
    #     '--where-field', dest='where_field', nargs='+', metavar='<field>=<regex>', required=False,
    #     help='Additional matches to header fields'
    # )
    # common.add_argument(
    #     '--out-fasta', dest='out_fasta', metavar='<filename>', required=False,
    #     help='A FASTA file (else writes to stdout)'
    # )
    # common.add_argument(
    #     '--out-metadata', dest='out_metadata', metavar='<filename>', required=False,
    #     help='A metadata file'
    # )
    # common.add_argument(
    #     '--header-append',  dest='header_append', nargs='*', metavar='<field>', required=False,
    #     help='Fields in the metadata table to construct the header'
    # )
    # common.add_argument(
    #     '--header-replace', dest='header_replace', nargs='*', metavar='<field>', required=False,
    #     help='Fields in the metadata table to construct the header'
    # )

    common.add_argument(
        "-v", "--verbose", dest="verbose", action="store_true",
        help="Run with high verbosity " "(debug level logging)",
    )
    common.add_argument(
        "--log-file", dest="log_file", metavar='<filename>', required=False, default=None,
        help="Log file to use (otherwise uses stdout, or stderr if out-fasta to stdout)"
    )
    # _______________________________  consensus  __________________________________#

    subparser_consensus = subparsers.add_parser(
        "consensus",
        parents=[common],
        help="Collapses sequences into consensus sequences based on grouping by index "
             "column or index field",
    )

    subparser_consensus.add_argument(
        '--in-fasta', dest='in_fasta', metavar='<filename>', required=False,
        help='One FASTA files of sequences (else reads from stdin)'
    )
    subparser_consensus.add_argument(
        '--in-metadata', dest='in_metadata', metavar='<filename>', required=True,
        help='CSV of metadata with same naming convention as fasta file'
    )
    subparser_consensus.add_argument(
        '--index-field', dest='index_field', metavar='<field>', required=True,
        help='Field(s) in the fasta header to match the metadata (else matches column names)'
    )
    subparser_consensus.add_argument(
        '--index-column', dest='index_column',default='header', metavar='<column>', required=False,
        help='Column(s) in the metadata file to use to match to the sequence'
    )
    subparser_consensus.add_argument(
        "--lineage", nargs='+', required=False, default="",
        help="Specific list of lineages to split by with others collpasing to nearest lineage."
    )
    subparser_consensus.add_argument(
        '--out-fasta', dest='out_fasta', metavar='<filename>', required=False, default="",
        help='A FASTA file (else writes to stdout)'
    )

    subparser_consensus.set_defaults(func=fastafunk.subcommands.consensus.run)

    # _______________________________  extract  __________________________________#

    subparser_extract = subparsers.add_parser(
        "extract",
        parents=[common],
        help="Extracts sequences based on matches to the metadata/tree",
    )

    subparser_extract.add_argument(
        '--in-fasta', dest='in_fasta', nargs='+', metavar='<filename>', required=True,
        help='One or more FASTA files of sequences (else reads from stdin)'
    )
    subparser_extract.add_argument(
        '--in-metadata', dest='in_metadata', nargs='+', metavar='<filename>', required=False,
        help='One or more CSV or TSV tables of metadata'
    )
    subparser_extract.add_argument(
        '--in-tree', dest='in_tree', nargs='+', metavar='<filename>', required=False,
        help='One or more tree files in either nexus or newick format'
    )
    subparser_extract.add_argument(
        '--out-fasta', dest='out_fasta', metavar='<filename>', required=False, default="",
        help='A FASTA file (else writes to stdout)'
    )
    subparser_extract.add_argument(
        '--reject-fasta', dest='reject_fasta', metavar='<filename>', required=False, default="",
        help='A FASTA file to write the omitted sequences '
    )
    subparser_extract.add_argument(
        '--low-memory', dest='low_memory', action='store_true', required=False,
        help='Extracts tip labels from trees using text wrangling instead of dendropy'
    )

    subparser_extract.set_defaults(func=fastafunk.subcommands.extract.run)

    # _______________________________  get_tips  __________________________________#

    subparser_get_tips = subparsers.add_parser(
        "get_tips",
        parents=[common],
        help="get_tips based on matches between the metadata/tree",
    )

    subparser_get_tips.add_argument(
        '--in-metadata', dest='in_metadata', nargs='+', metavar='<filename>', required=True,
        help='One or more CSV or TSV tables of metadata'
    )
    subparser_get_tips.add_argument(
        '--in-tree', dest='in_tree', nargs='+', metavar='<filename>', required=True,
        help='One or more tree files in either nexus or newick format'
    )
    subparser_get_tips.add_argument(
        '--out-tips', dest='out_tips', metavar='<filename>', required=True,
        help='A file to write the tip names '
    )
    subparser_get_tips.add_argument(
        '--low-memory', dest='low_memory', action='store_true', required=False,
        help='Extracts tip labels from trees using text wrangling instead of dendropy'
    )

    subparser_get_tips.set_defaults(func=fastafunk.subcommands.get_tips.run)

    # _______________________________  merge  __________________________________#

    subparser_merge = subparsers.add_parser(
        "merge",
        parents=[common],
        help="Merges two fasta files avoiding duplicates based on matches to "
             "metadata (takes the one in the first file)",
    )

    subparser_merge.add_argument(
        '--in-fasta', dest='in_fasta', nargs='+', metavar='<filename>', required=False,
        help='One or more FASTA files of sequences (else reads from stdin)'
    )
    subparser_merge.add_argument(
        '--in-metadata', dest='in_metadata', nargs='+', metavar='<filename>', required=True,
        help='One or more CSV or TSV tables of metadata'
    )
    subparser_merge.add_argument(
        '--index-column', dest='index_column', metavar='<column>', required=True,
        help='Column in the metadata file to use to match to the sequence'
    )
    subparser_merge.add_argument(
        '--out-metadata', dest='out_metadata', metavar='<filename>', required=True,
        help='A CSV file (else writes to stdout)'
    )
    subparser_merge.add_argument(
        '--out-fasta', dest='out_fasta', metavar='<filename>', required=False, default="",
        help='A FASTA file (else writes to stdout)'
    )
    subparser_merge.add_argument(
        '--low-memory', dest='low_memory', action='store_true', required=False,
        help='Assumes no duplicate sequences within a  FASTA so can use SeqIO index'
    )

    subparser_merge.set_defaults(func=fastafunk.subcommands.merge.run)

    # _______________________________  remove  __________________________________#

    subparser_remove = subparsers.add_parser(
        "remove",
        parents=[common],
        help="Removes sequences based on matches to the metadata",
    )

    subparser_remove.add_argument(
        '--in-fasta', dest='in_fasta', nargs='+', metavar='<filename>', required=True,
        help='One or more FASTA files of sequences (else reads from stdin)'
    )
    subparser_remove.add_argument(
        '--in-metadata', dest='in_metadata', nargs='+', metavar='<filename>', required=True,
        help='One or more CSV or TSV tables of metadata'
    )
    subparser_remove.add_argument(
        '--out-fasta', dest='out_fasta', metavar='<filename>', required=False, default="",
        help='A FASTA file (else writes to stdout)'
    )

    subparser_remove.set_defaults(func=fastafunk.subcommands.remove.run)

    # _______________________________  split  __________________________________#

    subparser_split = subparsers.add_parser(
        "split",
        parents=[common],
        help="Splits sequences into fasta files based on index column or index "
             "field",
    )

    subparser_split.add_argument(
        '--in-fasta', dest='in_fasta', metavar='<filename>', required=True,
        help='One FASTA files of sequences (else reads from stdin)'
    )
    subparser_split.add_argument(
        '--in-metadata', dest='in_metadata', metavar='<filename>', required=True,
        help='One CSV of metadata'
    )
    subparser_split.add_argument(
        '--index-column', dest='index_column', default='header', metavar='<column>', required=False,
        help='Column(s) in the metadata file to use to match to the sequence'
    )
    subparser_split.add_argument(
        '--index-field', dest='index_field', metavar='<field>', required=True,
        help='Field(s) in the fasta header to match the metadata (else matches column names)'
    )
    subparser_split.add_argument(
        "--lineage", nargs='+', required=False, default="",
        help="Specific list of lineages to split by with others collpasing to nearest lineage."
    )
    subparser_split.add_argument(
        "--lineage-csv", dest='lineage_csv', required=False, default="",
        help="CSV with lineage and outgroup columns defining the lineages to split by."
    )
    subparser_split.add_argument(
        "--aliases", dest='aliases', required=False, default="",
        help="JSON with aliases for lettered lineages."
    )
    subparser_split.add_argument(
        '--out-folder', dest='out_prefix', metavar='<filename>', default="./", required=False,
        help='A directory for output FASTA files'
    )

    subparser_split.set_defaults(func=fastafunk.subcommands.split.run)

    # _______________________________  count  __________________________________#

    subparser_count = subparsers.add_parser(
        "count",
        parents=[common],
        help="Counts the number in each group, defined by a number of metadata columns",
    )

    subparser_count.add_argument(
        '--in-metadata', dest='in_metadata', nargs='+', metavar='<filename>', required=True,
        help='One or more CSV or TSV tables of metadata'
    )
    subparser_count.add_argument(
        '--group-column', dest='group_column', nargs='+', metavar='<column>', required=True,
        help='Column(s) in the metadata file to define groupings'
    )

    subparser_count.set_defaults(func=fastafunk.subcommands.count.run)

    # _______________________________  subsample  __________________________________#

    subparser_subsample = subparsers.add_parser(
        "subsample",
        parents=[common],
        help="Subsamples a fasta based on groups defined by a metadata file",
    )

    subparser_subsample.add_argument(
        '--in-fasta', dest='in_fasta', nargs='+', metavar='<filename>', required=False,
        help='One or more FASTA files of sequences (else reads from stdin)'
    )
    subparser_subsample.add_argument(
        '--in-metadata', dest='in_metadata', nargs='+', metavar='<filename>', required=True,
        help='One or more CSV or TSV tables of metadata'
    )
    subparser_subsample.add_argument(
        '--index-field', dest='index_field', nargs='+', metavar='<field>', required=False,
        help='Field(s) in the fasta header to match the metadata (else matches column names)'
    )
    subparser_subsample.add_argument(
        '--index-column', dest='index_column', metavar='<column>', required=True,
        help='Column in the metadata file to use to match to the sequence'
    )
    subparser_subsample.add_argument(
        '--group-column', dest='group_column', nargs='+', metavar='<column>', required=True,
        help='Column(s) in the metadata file to define groupings'
    )
    subparser_subsample.add_argument(
        '--where-field', dest='where_field', metavar='<field>=<regex>', required=False,
        help='Additional matches to header fields'
    )
    subparser_subsample.add_argument(
        '--out-fasta', dest='out_fasta', metavar='<filename>', required=False, default="",
        help='A FASTA file (else writes to stdout)'
    )
    subparser_subsample.add_argument(
        '--out-metadata', dest='out_metadata', metavar='<filename>', required=False,
        help='A metadata file'
    )
    subparser_subsample.add_argument(
        "--target-file", dest='target_file', metavar="<filename>", required=False, default="",
        help="CSV file of target numbers per group e.g. an edited version of the count output"
    )
    subparser_subsample.add_argument(
        '--select-by-max-column', dest='select_by_max_column', metavar='<column>', required=False,
        help='Column in the metadata file maximize over when subsetting'
    )
    subparser_subsample.add_argument(
        '--select-by-min-column', dest='select_by_min_column', metavar='<column>', required=False,
        help='Column in the metadata file minimize over when subsetting'
    )
    subparser_subsample.add_argument(
        "--sample-size", dest='sample_size', metavar="<int>", type=int, required=False, default=10,
        help="The number of samples per group to select if not specified by target file"
    )
    subparser_subsample.add_argument(
        "--exclude-uk", dest='exclude_uk', action="store_true",
        help="Includes all UK samples in subsample, and additionally keeps the target number of "
             "non-UK samples per group"
    )
    subparser_subsample.add_argument(
        '--low-memory', dest='low_memory', action='store_true', required=False,
        help='Assumes no duplicate sequences within a  FASTA so can use SeqIO index'
    )

    subparser_subsample.set_defaults(func=fastafunk.subcommands.subsample.run)

    # _______________________________  annotate  __________________________________#

    subparser_annotate = subparsers.add_parser(
        "annotate",
        parents=[common],
        help="Subsamples a fasta based on groups defined by a metadata file",
    )

    subparser_annotate.add_argument(
        '--in-fasta', dest='in_fasta', nargs='+', metavar='<filename>', required=False,
        help='One or more FASTA files of sequences (else reads from stdin)'
    )
    subparser_annotate.add_argument(
        '--in-metadata', dest='in_metadata', nargs='+', metavar='<filename>', required=False,
        help='One or more CSV or TSV tables of metadata'
    )
    subparser_annotate.add_argument(
        '--index-field', dest='index_field', nargs='+', metavar='<field>', required=False,
        help='Field(s) in the fasta header to match the metadata (else matches column names)'
    )
    subparser_annotate.add_argument(
        '--index-column', dest='index_column', metavar='<column>', required=True,
        help='Column in the metadata file to use to match to the sequence'
    )
    subparser_annotate.add_argument(
        '--out-fasta', dest='out_fasta', metavar='<filename>', required=False, default="",
        help='A FASTA file (else writes to stdout)'
    )
    subparser_annotate.add_argument(
        '--out-metadata', dest='out_metadata', metavar='<filename>', required=False,
        help='A metadata file'
    )
    subparser_annotate.add_argument(
        '--header-delimiter', dest='header_delimiter', default='|', metavar='<symbol>', required=False,
        help='Header delimiter'
    )
    subparser_annotate.add_argument(
        '--add-cov-id', dest='add_cov_id', action="store_true", required=False,
        help='Parses header for COG or GISAID unique id and stores'
    )
    subparser_annotate.add_argument(
        '--low-memory', dest='low_memory', action='store_true', required=False,
        help='Assumes no duplicate sequences within a  FASTA so can use SeqIO index'
    )

    subparser_annotate.set_defaults(func=fastafunk.subcommands.annotate.run)

    # _______________________________  unwrap  __________________________________#

    subparser_unwrap = subparsers.add_parser(
        "unwrap",
        parents=[common],
        help="Remove whitespace from the sequences in 2 line fasta format"
    )

    subparser_unwrap.add_argument(
        '--in-fasta', dest='in_fasta', nargs='+', metavar='<filename>', required=True,
        help='One or more FASTA files of sequences (else reads from stdin)'
    )
    subparser_unwrap.add_argument(
        '--out-fasta', dest='out_fasta', metavar='<filename>', required=False, default="",
        help='A FASTA file (else writes to stdout)'
    )

    subparser_unwrap.set_defaults(func=fastafunk.subcommands.unwrap.run)


    # _______________________________  strip  __________________________________#

    subparser_strip = subparsers.add_parser(
        "strip",
        parents=[common],
        help="Strip sites based on various options"
    )

    subparser_strip.add_argument(
        '--in-fasta', dest='in_fasta', nargs='+', metavar='<filename>', required=True,
        help='One or more FASTA files of sequences (else reads from stdin)'
    )
    subparser_strip.add_argument(
        '--gap', dest='gap', default=False, action='store_true', required=False,
        help='Remove gaps from sequences (Default:False)'
    )
    subparser_strip.add_argument(
        '--ambiguity', dest='ambiguity', default=False, action='store_true', required=False,
        help='Remove ambiguous sites from sequences ("N") (Default:False)'
    )
    subparser_strip.add_argument(
        '--missing', dest='missing', default=False, action='store_true', required=False,
        help='Remove missing sites from sequences ("?") (Default:False)'
    )
    subparser_strip.add_argument(
        '--keep-alignment', dest='keep_alignment', default=False, action='store_true', required=False,
        help='Remove gaps shared by all sequences at the same site (Default:False)'
    )
    subparser_strip.add_argument(
        '--front', dest='front', default=False, action='store_true', required=False,
        help='Remove only from the front of the sequence (Default:False)'
    )
    subparser_strip.add_argument(
        '--back', dest='back', default=False, action='store_true', required=False,
        help='Remove only from the back of the sequence (Default:False)'
    )
    subparser_strip.add_argument(
        '--out-fasta', dest='out_fasta', metavar='<filename>', required=False, default="stripped.fasta",
        help='A FASTA file (else writes to stdout)'
    )

    subparser_strip.set_defaults(func=fastafunk.subcommands.strip.run)


    # __________________________________________________________________________#
    # _______________________________  add_columns  __________________________________#

    subparser_add_columns = subparsers.add_parser(
        "add_columns",
        parents=[common],
        help="Add columns to a metadata file from row matches in another file"
    )

    subparser_add_columns.add_argument(
        '--in-metadata', dest='in_metadata', metavar='<filename>', required=True,
        help='ONE CSV table of metadata'
    )

    subparser_add_columns.add_argument(
        '--in-data', dest='in_data', metavar='<filename>', required=True,
        help='One CSV table of additional data. Must have --index-column in common with metadata'
    )

    subparser_add_columns.add_argument(
        '--index-column', dest='index_column', metavar='<column>', required=True,
        help='Column in the metadata files used to match rows'
    )

    subparser_add_columns.add_argument(
        '--join-on', dest='join_on', metavar='<column>', required=True,
        help='Column in the data file used to match rows'
    )

    subparser_add_columns.add_argument(
        '--new-columns', dest='new_columns', nargs='+', metavar='<column>', required=False,
        help='Column(s) in the in_data file to add to the metadata, if not provided, all columns added'
    )

    subparser_add_columns.add_argument(
        '--out-metadata', dest='out_metadata', metavar='<filename>', required=True,
        help='A metadata file to write'
    )
    subparser_add_columns.add_argument(
        '--where-column', dest='where_column', nargs='+', metavar='<column>=<regex>', required=False,
        help='Additional matches to columns e.g. if want to rename'
    )
    subparser_add_columns.add_argument('--force-overwrite', dest='force_overwrite', action='store_true', required=False,
                                       help='Overwrite even if new data is blank/None')

    subparser_add_columns.set_defaults(func=fastafunk.subcommands.add_columns.run)

    # _______________________________  drop_columns  __________________________________#

    subparser_drop_columns = subparsers.add_parser(
        "drop_columns",
        parents=[common],
        help="Drop columns from a metadata file "
    )

    subparser_drop_columns.add_argument(
        '--in-metadata', dest='in_metadata', metavar='<filename>', required=True,
        help='ONE CSV table of metadata'
    )

    subparser_drop_columns.add_argument(
        '--columns', dest='columns', nargs='+', metavar='<column>', required=False,
        help='Column(s) in the metadata to drop'
    )

    subparser_drop_columns.add_argument(
        '--out-metadata', dest='out_metadata', metavar='<filename>', required=True,
        help='A metadata file to write'
    )

    subparser_drop_columns.set_defaults(func=fastafunk.subcommands.drop_columns.run)

    # _______________________________  filter_column  __________________________________#

    subparser_filter_column = subparsers.add_parser(
        "filter_column",
        parents=[common],
        help="Filter metadata file based on column"
    )

    subparser_filter_column.add_argument(
        '--in-metadata', dest='in_metadata', metavar='<filename>', required=True,
        help='ONE table of metadata'
    )

    subparser_filter_column.add_argument(
        '--column', dest='column', metavar='<column>', required=True,
        help='Column in the metadata to filter on'
    )

    subparser_filter_column.add_argument(
        '--out-metadata', dest='out_metadata', metavar='<filename>', required=True,
        help='A metadata file to write'
    )

    subparser_filter_column.add_argument('--is_true', action='store_true', required=False, help='filter if column is true')
    subparser_filter_column.add_argument('--is_false', action='store_true', required=False, help='filter if column is false')

    subparser_filter_column.set_defaults(func=fastafunk.subcommands.filter_column.run)


    # _______________________________  new  __________________________________#

    subparser_new = subparsers.add_parser(
        "new",
        parents=[common],
        help="Identifies new or updated entries in metadata file and associated fasta entries",
    )

    subparser_new.add_argument(
        '--in-fasta', dest='in_fasta', nargs='+', metavar='<filename>', required=False,
        help='One or more FASTA files of sequences (else reads from stdin)'
    )
    subparser_new.add_argument(
        '--in-metadata', dest='in_metadata', nargs='+', metavar='<filename>', required=True,
        help='One or more CSV or TSV tables of metadata'
    )
    subparser_new.add_argument(
        '--index-column', dest='index_column', metavar='<column>', required=True,
        help='Column in the metadata file to use to match to the sequence'
    )
    subparser_new.add_argument(
        '--date-column', dest='date_column', metavar='<column>', required=True,
        help='Column in the metadata file containing date values to use for filtering'
    )
    subparser_new.add_argument(
        '--out-fasta', dest='out_fasta', metavar='<filename>', required=False, default="",
        help='A FASTA file (else writes to stdout)'
    )
    subparser_new.add_argument(
        '--out-metadata', dest='out_metadata', metavar='<filename>', required=False,
        help='A metadata file'
    )

    subparser_new.set_defaults(func=fastafunk.subcommands.new.run)

    # _______________________________  fetch  __________________________________#

    subparser_fetch = subparsers.add_parser(
        "fetch",
        parents=[common],
        help="Fetches sequences within fasta files which have a corresponding metadata entry, "
             "keeping the last where there are duplicates",
    )

    subparser_fetch.add_argument(
        '--in-fasta', dest='in_fasta', nargs='+', metavar='<filename>', required=False,
        help='One or more FASTA files of sequences (else reads from stdin)'
    )
    subparser_fetch.add_argument(
        '--in-metadata', dest='in_metadata', metavar='<filename>', required=True,
        help='CSV or TSV of metadata with same naming convention as fasta file'
    )
    subparser_fetch.add_argument(
        '--index-column', dest='index_column', metavar='<column>', required=True,
        help='Column in the metadata file to use to match to the sequence'
    )
    subparser_fetch.add_argument(
        '--filter-column', dest='filter_column', nargs='+', metavar='<column>', required=False,
        help='Metadata column name(s) to keep'
    )
    subparser_fetch.add_argument(
        '--where-column', dest='where_column', nargs='+', metavar='<column>=<regex>', required=False,
        help='Additional matches to columns e.g. if want to rename'
    )
    subparser_fetch.add_argument(
        '--restrict', dest='restrict', action='store_true', required=False,
        help='Only outputs metadata rows with a corresponding fasta entry'
    )
    subparser_fetch.add_argument(
        '--keep-omit-rows', dest='keep_omit_rows', action='store_true', required=False,
        help='Allows rows with with metadata saying omit to be kept'
    )
    subparser_fetch.add_argument(
        '--low-memory', dest='low_memory', action='store_true', required=False,
        help='Assumes no duplicate sequences within a  FASTA so can use SeqIO index'
    )
    subparser_fetch.add_argument(
        '--out-fasta', dest='out_fasta', metavar='<filename>', required=False, default="",
        help='A FASTA file (else writes to stdout)'
    )
    subparser_fetch.add_argument(
        '--out-metadata', dest='out_metadata', metavar='<filename>', required=False,
        help='A metadata file'
    )

    subparser_fetch.set_defaults(func=fastafunk.subcommands.fetch.run)

    # _______________________________  shuffle  __________________________________#

    subparser_shuffle = subparsers.add_parser(
        "shuffle",
        parents=[common],
        help="Shuffles lines of a metadata file",
    )

    subparser_shuffle.add_argument(
        '--in-metadata', dest='in_metadata', metavar='<filename>', required=True,
        help='CSV or TSV of metadata'
    )

    subparser_shuffle.add_argument(
        '--out-metadata', dest='out_metadata', metavar='<filename>', required=True,
        help='CSV or TSV of metadata'
    )

    subparser_shuffle.set_defaults(func=fastafunk.subcommands.shuffle.run)

    # _________________________________________________________________________#


    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
