"""
This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import unittest
import filecmp
import glob

from fastafunk.fetch import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'fetch')

class TestFetch(unittest.TestCase):
    def test_run_fetch(self):
        in_fasta = ["%s/in_fasta.fa" % data_dir]
        in_metadata = "%s/in_metadata.csv" %data_dir
        index_column = "first"
        filter_column = ["header", "first", "id", "BLAH", "parrot", "why_excluded"]
        where_column = None
        restrict = False
        out_fasta = "%s/tmp.fetch.fa" % data_dir
        out_metadata = "%s/tmp.fetch.csv" % data_dir
        log_file = "%s/tmp.fetch.log" % data_dir
        low_memory = False
        keep_omit_rows = True
        fetch_fasta(in_fasta, in_metadata, index_column, filter_column, where_column, restrict, out_fasta,out_metadata, log_file, low_memory, keep_omit_rows)
        expected_metadata = "%s/expected_fetch.csv" %data_dir
        expected_fasta = "%s/expected_fetch.fa" % data_dir
        self.assertTrue(filecmp.cmp(out_metadata, expected_metadata, shallow=False))
        self.assertTrue(filecmp.cmp(out_fasta, expected_fasta, shallow=False))
        os.unlink(out_metadata)
        os.unlink(out_fasta)
        os.unlink(log_file)

    def test_run_fetch_reorder(self):
        in_fasta = ["%s/in_fasta.fa" % data_dir]
        in_metadata = "%s/in_metadata.csv" %data_dir
        index_column = "first"
        filter_column = ["header", "first", "id", "BLAH", "parrot", "why_excluded", "edin_omit"]
        where_column = None
        restrict = False
        out_fasta = "%s/tmp.fetch_reorder.fa" % data_dir
        out_metadata = "%s/tmp.fetch_reorder.csv" % data_dir
        log_file = "%s/tmp.fetch_reorder.log" % data_dir
        low_memory = False
        keep_omit_rows = True
        fetch_fasta(in_fasta, in_metadata, index_column, filter_column, where_column, restrict, out_fasta,out_metadata, log_file, low_memory, keep_omit_rows)
        expected_metadata = "%s/expected_fetch_reorder.csv" %data_dir
        expected_fasta = "%s/expected_fetch.fa" % data_dir
        self.assertTrue(filecmp.cmp(out_metadata, expected_metadata, shallow=False))
        self.assertTrue(filecmp.cmp(out_fasta, expected_fasta, shallow=False))
        os.unlink(out_metadata)
        os.unlink(out_fasta)
        os.unlink(log_file)

    def test_run_fetch_restrict(self):
        in_fasta = ["%s/in_fasta.fa" % data_dir]
        in_metadata = "%s/in_metadata.csv" %data_dir
        index_column = "first"
        filter_column = ["header", "first", "id", "BLAH", "parrot", "why_excluded", "edin_omit"]
        where_column = None
        restrict = True
        out_fasta = "%s/tmp.fetch_restrict.fa" % data_dir
        out_metadata = "%s/tmp.fetch_restrict.csv" % data_dir
        log_file = "%s/tmp.fetch_restrict.log" % data_dir
        low_memory = False
        keep_omit_rows = True
        fetch_fasta(in_fasta, in_metadata, index_column, filter_column, where_column, restrict, out_fasta,out_metadata, log_file, low_memory, keep_omit_rows)
        expected_metadata = "%s/expected_fetch_restrict.csv" %data_dir
        expected_fasta = "%s/expected_fetch.fa" % data_dir
        self.assertTrue(filecmp.cmp(out_metadata, expected_metadata, shallow=False))
        self.assertTrue(filecmp.cmp(out_fasta, expected_fasta, shallow=False))
        os.unlink(out_metadata)
        os.unlink(out_fasta)
        os.unlink(log_file)

    def test_run_fetch_low_memory(self):
        in_fasta = ["%s/in_fasta.fa" % data_dir]
        in_metadata = "%s/in_metadata.csv" %data_dir
        index_column = "first"
        filter_column = ["header", "first", "id", "BLAH", "parrot", "why_excluded", "edin_omit"]
        where_column = None
        restrict = True
        out_fasta = "%s/tmp.fetch_low_memory.fa" % data_dir
        out_metadata = "%s/tmp.fetch_low_memory.csv" % data_dir
        log_file = "%s/tmp.fetch_low_memory.log" % data_dir
        low_memory = True
        keep_omit_rows = True
        fetch_fasta(in_fasta, in_metadata, index_column, filter_column, where_column, restrict, out_fasta,out_metadata, log_file, low_memory, keep_omit_rows)
        expected_metadata = "%s/expected_fetch_restrict.csv" %data_dir
        expected_fasta = "%s/expected_fetch.fa" % data_dir
        self.assertTrue(filecmp.cmp(out_metadata, expected_metadata, shallow=False))
        self.assertTrue(filecmp.cmp(out_fasta, expected_fasta, shallow=False))
        os.unlink(out_metadata)
        os.unlink(out_fasta)
        os.unlink(log_file)

    def test_run_fetch_omit_rows(self):
        in_fasta = ["%s/in_fasta.fa" % data_dir]
        in_metadata = "%s/in_metadata.csv" %data_dir
        index_column = "first"
        filter_column = ["header", "first", "id", "BLAH", "parrot", "why_excluded", "edin_omit"]
        where_column = None
        restrict = True
        out_fasta = "%s/tmp.fetch_omit_rows.fa" % data_dir
        out_metadata = "%s/tmp.fetch_omit_rows.csv" % data_dir
        log_file = "%s/tmp.fetch_omit_rows.log" % data_dir
        low_memory = True
        keep_omit_rows = False
        fetch_fasta(in_fasta, in_metadata, index_column, filter_column, where_column, restrict, out_fasta,out_metadata, log_file, low_memory, keep_omit_rows)
        expected_metadata = "%s/expected_fetch_omit_rows.csv" %data_dir
        expected_fasta = "%s/expected_fetch_omit_rows.fa" % data_dir
        self.assertTrue(filecmp.cmp(out_metadata, expected_metadata, shallow=False))
        self.assertTrue(filecmp.cmp(out_fasta, expected_fasta, shallow=False))
        os.unlink(out_metadata)
        os.unlink(out_fasta)
        os.unlink(log_file)

    def test_run_fetch_omit_rows_no_restrict(self):
        in_fasta = ["%s/in_fasta.fa" % data_dir]
        in_metadata = "%s/in_metadata.csv" %data_dir
        index_column = "first"
        filter_column = ["header", "first", "id", "BLAH", "parrot", "why_excluded", "edin_omit"]
        where_column = None
        restrict = False
        out_fasta = "%s/tmp.fetch_omit_rows_no_restrict.fa" % data_dir
        out_metadata = "%s/tmp.fetch_omit_rows_no_restrict.csv" % data_dir
        log_file = "%s/tmp.fetch_omit_rows_no_restrict.log" % data_dir
        low_memory = True
        keep_omit_rows = False
        fetch_fasta(in_fasta, in_metadata, index_column, filter_column, where_column, restrict, out_fasta,out_metadata, log_file, low_memory, keep_omit_rows)
        expected_metadata = "%s/expected_fetch_omit_rows_no_restrict.csv" %data_dir
        expected_fasta = "%s/expected_fetch_omit_rows.fa" % data_dir
        self.assertTrue(filecmp.cmp(out_metadata, expected_metadata, shallow=False))
        self.assertTrue(filecmp.cmp(out_fasta, expected_fasta, shallow=False))
        os.unlink(out_metadata)
        os.unlink(out_fasta)
        os.unlink(log_file)

    def test_run_fetch_no_filter_column(self):
        in_fasta = ["%s/in_fasta.fa" % data_dir]
        in_metadata = "%s/in_metadata.csv" %data_dir
        index_column = "first"
        filter_column = None
        where_column = None
        restrict = False
        out_fasta = "%s/tmp.fetch_no_filter_column.fa" % data_dir
        out_metadata = "%s/tmp.fetch_no_filter_column.csv" % data_dir
        log_file = "%s/tmp.fetch_no_filter_column.log" % data_dir
        low_memory = True
        keep_omit_rows = True
        fetch_fasta(in_fasta, in_metadata, index_column, filter_column, where_column, restrict, out_fasta,out_metadata, log_file, low_memory, keep_omit_rows)
        expected_metadata = "%s/in_metadata.csv" %data_dir
        expected_fasta = "%s/expected_fetch.fa" % data_dir
        self.assertTrue(filecmp.cmp(out_metadata, expected_metadata, shallow=False))
        self.assertTrue(filecmp.cmp(out_fasta, expected_fasta, shallow=False))
        os.unlink(out_metadata)
        os.unlink(out_fasta)
        os.unlink(log_file)
