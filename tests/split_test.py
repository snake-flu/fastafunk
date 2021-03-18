"""
This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import unittest
import filecmp
import glob

from fastafunk.split import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'split')

class TestSplit(unittest.TestCase):
    def test_run_split(self):
        in_fasta = "%s/in_fasta.fa" % data_dir
        in_metadata = "%s/in_metadata.csv" %data_dir
        index_field = "lineage"
        index_column = "sequence_name"
        lineage = None
        lineage_csv = "%s/lineage_splits.csv" %data_dir
        aliases = "%s/lineage_aliases.json" %data_dir
        out_folder = "%s/tmp.split_" %data_dir
        log_file = "%s/tmp.split.log" % data_dir
        split_fasta(in_fasta, in_metadata, index_field, index_column, lineage, lineage_csv, aliases, out_folder, log_file)
        for lineage in ["A", "B", "B.1", "B.1.1.7"]:
            expected_fasta = "%s/expected_%s.fasta" % (data_dir, lineage)
            out_fasta = "%s/tmp.split_%s.fasta" % (data_dir, lineage)
            self.assertTrue(filecmp.cmp(out_fasta, expected_fasta, shallow=False))
            os.unlink(out_fasta)
        os.unlink(log_file)

