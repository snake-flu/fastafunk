"""
This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import unittest
import filecmp
import glob

from fastafunk.add_columns import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'add_columns')

class TestAddColumns(unittest.TestCase):
    def test_run_add_columns(self):
        in_metadata = "%s/simple_metadata.csv" %data_dir
        in_data = "%s/new_metadata.csv" %data_dir
        index_column = "header"
        join_on = "taxon"
        new_columns = ["mut1","parrot"]
        out_metadata = "%s/tmp.add_columns.csv" %data_dir
        where_column = None
        log_file = "%s/tmp.add_columns.log" %data_dir
        expected = "%s/expected_add_columns.csv" %data_dir
        add_columns(in_metadata, in_data, index_column, join_on, new_columns, out_metadata, where_column, log_file)
        self.assertTrue(filecmp.cmp(out_metadata, expected, shallow=False))
        os.unlink(out_metadata)
        os.unlink(log_file)

    def test_run_add_columns_no_new_columns(self):
        in_metadata = "%s/simple_metadata.csv" %data_dir
        in_data = "%s/new_metadata.csv" %data_dir
        index_column = "header"
        join_on = "taxon"
        new_columns = None
        out_metadata = "%s/tmp.add_columns1.csv" %data_dir
        where_column = None
        log_file = "%s/tmp.add_columns.log" %data_dir
        expected = "%s/expected_add_columns2.csv" %data_dir
        add_columns(in_metadata, in_data, index_column, join_on, new_columns, out_metadata, where_column, log_file)
        self.assertTrue(filecmp.cmp(out_metadata, expected, shallow=False))
        os.unlink(out_metadata)
        os.unlink(log_file)

    def test_run_add_columns_no_new_columns2(self):
        in_metadata = "%s/simple_metadata.csv" %data_dir
        in_data = "%s/new_metadata.csv" %data_dir
        index_column = "header"
        join_on = "taxon"
        new_columns = []
        out_metadata = "%s/tmp.add_columns2.csv" %data_dir
        where_column = None
        log_file = "%s/tmp.add_columns.log" %data_dir
        expected = "%s/expected_add_columns2.csv" %data_dir
        add_columns(in_metadata, in_data, index_column, join_on, new_columns, out_metadata, where_column, log_file)
        self.assertTrue(filecmp.cmp(out_metadata, expected, shallow=False))
        os.unlink(out_metadata)
        os.unlink(log_file)

    def test_run_add_columns_force_overwrite(self):
        in_metadata = "%s/simple_metadata.csv" %data_dir
        in_data = "%s/new_metadata.csv" %data_dir
        index_column = "header"
        join_on = "taxon"
        new_columns = ["mut1","parrot"]
        out_metadata = "%s/tmp.add_columns.csv" %data_dir
        where_column = None
        log_file = "%s/tmp.add_columns.log" %data_dir
        force_overwrite = True
        expected = "%s/expected_add_columns3.csv" %data_dir
        add_columns(in_metadata, in_data, index_column, join_on, new_columns, out_metadata, where_column, log_file, force_overwrite)
        self.assertTrue(filecmp.cmp(out_metadata, expected, shallow=False))
        os.unlink(out_metadata)
        os.unlink(log_file)
