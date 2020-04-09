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
    # common.add_argument(
    #     '--header-delimiter', dest='header_delimiter', default='|', metavar='<symbol>', required=False,
    #     help='Header delimiter'
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
        '--in-fasta', dest='in_fasta', nargs='+', metavar='<filename>', required=False,
        help='One or more FASTA files of sequences (else reads from stdin)'
    )
    subparser_consensus.add_argument(
        '--in-metadata', dest='in_metadata', metavar='<filename>', required=True,
        help='CSV of metadata with same naming convention as fasta file'
    )
    subparser_consensus.add_argument(
        '--index-field', dest='index_field', nargs='+', metavar='<field>', required=True,
        help='Field(s) in the fasta header to match the metadata (else matches column names)'
    )
    subparser_consensus.add_argument(
        '--index-column', dest='index_column', nargs='+', metavar='<column>', required=False,
        help='Column(s) in the metadata file to use to match to the sequence'
    )
    subparser_consensus.add_argument(
        '--out-fasta', dest='out_fasta', metavar='<filename>', required=False, default="",
        help='A FASTA file (else writes to stdout)'
    )
    subparser_consensus.add_argument(
        "--clade-file", required=False, default="", type=str, nargs='+', metavar="<filename>",
        help="Text file including specific traits to collapse by"
    )

    subparser_consensus.set_defaults(func=fastafunk.subcommands.consensus.run)

    # _______________________________  extract  __________________________________#

    subparser_extract = subparsers.add_parser(
        "extract",
        parents=[common],
        help="Extracts sequences based on matches to the metadata",
    )

    subparser_extract.add_argument(
        '--in-fasta', dest='in_fasta', nargs='+', metavar='<filename>', required=True,
        help='One or more FASTA files of sequences (else reads from stdin)'
    )
    subparser_extract.add_argument(
        '--in-metadata', dest='in_metadata', nargs='+', metavar='<filename>', required=True,
        help='One or more CSV or TSV tables of metadata'
    )
    subparser_extract.add_argument(
        '--out-fasta', dest='out_fasta', metavar='<filename>', required=False, default="",
        help='A FASTA file (else writes to stdout)'
    )

    subparser_extract.set_defaults(func=fastafunk.subcommands.extract.run)

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
        '--out-fasta', dest='out_fasta', metavar='<filename>', required=False, default="",
        help='A FASTA file (else writes to stdout)'
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
        '--in-fasta', dest='in_fasta', nargs='+', metavar='<filename>', required=True,
        help='One or more FASTA files of sequences (else reads from stdin)'
    )
    subparser_split.add_argument(
        '--in-metadata', dest='in_metadata', metavar='<filename>', required=True,
        help='One CSV of metadata'
    )
    subparser_split.add_argument(
        '--index-column', dest='index_column', nargs='+', metavar='<column>', required=False,
        help='Column(s) in the metadata file to use to match to the sequence'
    )
    subparser_split.add_argument(
        '--index-field', dest='index_field', nargs='+', metavar='<field>', required=True,
        help='Field(s) in the fasta header to match the metadata (else matches column names)'
    )
    subparser_split.add_argument(
        '--out-folder', dest='out_folder', metavar='<filename>', default="./", required=False,
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
        '--index-column', dest='index_column', nargs='+', metavar='<column>', required=True,
        help='Column(s) in the metadata file to use to match to the sequence'
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
        '--index-column', dest='index_column', nargs='+', metavar='<column>', required=True,
        help='Column(s) in the metadata file to use to match to the sequence'
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
        "--sample-size", dest='sample_size', metavar="<int>", type=int, required=False, default=10,
        help="The number of samples per group to select if not specified by target file"
    )
    subparser_subsample.add_argument(
        "--exclude-uk", dest='exclude_uk', action="store_true",
        help="Includes all UK samples in subsample, and additionally keeps the target number of "
             "non-UK samples per group"
    )

    subparser_subsample.set_defaults(func=fastafunk.subcommands.subsample.run)
    # ___________________________________________________________________________#


    args = parser.parse_args()

    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()