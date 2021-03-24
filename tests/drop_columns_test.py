"""
This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import unittest
import filecmp
import glob

from fastafunk.drop_columns import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'drop_columns')

class TestAddColumns(unittest.TestCase):
    def test_run_drop_columns(self):
        in_metadata = "%s/simple_metadata.csv" %data_dir
        columns = ["id"]
        out_metadata = "%s/tmp.drop_columns.csv" %data_dir
        log_file = "%s/tmp.drop_columns.log" %data_dir
        expected = "%s/expected_drop_columns1.csv" %data_dir
        drop_columns(in_metadata, columns, out_metadata, log_file)
        self.assertTrue(filecmp.cmp(out_metadata, expected, shallow=False))
        os.unlink(out_metadata)
        os.unlink(log_file)

    def test_run_drop_columns2(self):
        in_metadata = "%s/simple_metadata.csv" %data_dir
        columns = ["header","first","BLAH"]
        out_metadata = "%s/tmp.drop_columns2.csv" %data_dir
        log_file = "%s/tmp.drop_columns2.log" %data_dir
        expected = "%s/expected_drop_columns2.csv" %data_dir
        drop_columns(in_metadata, columns, out_metadata, log_file)
        self.assertTrue(filecmp.cmp(out_metadata, expected, shallow=False))
        os.unlink(out_metadata)
        os.unlink(log_file)
