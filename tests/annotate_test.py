"""
This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import os
import unittest
import filecmp
import pandas as pd
import pandas.testing as pd_testing
from Bio import SeqIO

from fastafunk.annotate import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'annotate')

class TestAnnotate(unittest.TestCase):
    def assertDataframeEqual(self, a, b, msg):
        try:
            pd_testing.assert_frame_equal(a, b)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def setUp(self):
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)

    def test_annotate_simple(self):
        in_fasta = ["%s/annotate.fasta" %data_dir]
        in_metadata = ["%s/simple_metadata.csv" % data_dir]
        index_column = None
        index_field = None
        out_fasta = "%s/tmp.annotated.fasta" %data_dir
        out_metadata = None
        header_delimiter = "|"
        log_file = None
        annotate(in_fasta, in_metadata, index_column, index_field, out_fasta, out_metadata, header_delimiter, log_file)

        expected = "%s/expect_simple.fasta" %data_dir
        self.assertTrue(filecmp.cmp(out_fasta, expected, shallow=False))
        os.unlink(out_fasta)

    def test_annotate_simple_index_column_index_field(self):
        in_fasta = ["%s/annotate.fasta" %data_dir]
        in_metadata = ["%s/simple_metadata.csv" % data_dir]
        index_column = ""
        index_field = None
        out_fasta = "%s/tmp.annotated.fasta" %data_dir
        out_metadata = None
        header_delimiter = "|"
        log_file = None
        annotate(in_fasta, in_metadata, index_column, index_field, out_fasta, out_metadata, header_delimiter, log_file)

        expected = "%s/expect_simple.fasta" %data_dir
        self.assertTrue(filecmp.cmp(out_fasta, expected, shallow=False))
        os.unlink(out_fasta)

    def test_annotate_metadata_simple(self):
        in_fasta = ["%s/annotate.fasta" %data_dir]
        in_metadata = ["%s/simple_metadata.csv" % data_dir]
        index_column = ""
        index_field = None
        out_fasta = "%s/tmp.annotated.fasta" %data_dir
        out_metadata = "%s/tmp.annotated.csv" %data_dir
        header_delimiter = "|"
        log_file = None
        annotate(in_fasta, in_metadata, index_column, index_field, out_fasta, out_metadata, header_delimiter, log_file)

        expected = "%s/expect_simple.fasta" % data_dir
        self.assertTrue(filecmp.cmp(out_fasta, expected, shallow=False))
        os.unlink(out_fasta)

        expected = "%s/expect_annotated.csv" %data_dir
        self.assertTrue(filecmp.cmp(out_metadata, expected, shallow=False))
        os.unlink(out_metadata)

    def test_annotate_metadata_superset(self):
        in_fasta = ["%s/annotate.fasta" %data_dir]
        in_metadata = ["%s/superset_metadata.csv" % data_dir]
        index_column = ""
        index_field = None
        out_fasta = "%s/tmp.annotated.fasta" %data_dir
        out_metadata = "%s/tmp.annotated.csv" %data_dir
        header_delimiter = "|"
        log_file = None
        annotate(in_fasta, in_metadata, index_column, index_field, out_fasta, out_metadata, header_delimiter, log_file)

        expected = "%s/expect_simple.fasta" % data_dir
        self.assertTrue(filecmp.cmp(out_fasta, expected, shallow=False))
        os.unlink(out_fasta)

        expected = "%s/expect_annotated_superset.csv" %data_dir
        self.assertTrue(filecmp.cmp(out_metadata, expected, shallow=False))
        os.unlink(out_metadata)

    def test_annotate_metadata_subset(self):
        in_fasta = ["%s/annotate.fasta" %data_dir]
        in_metadata = ["%s/subset_metadata.csv" % data_dir]
        index_column = ""
        index_field = None
        out_fasta = "%s/tmp.annotated.fasta" %data_dir
        out_metadata = "%s/tmp.annotated.csv" %data_dir
        header_delimiter = "|"
        log_file = None
        annotate(in_fasta, in_metadata, index_column, index_field, out_fasta, out_metadata, header_delimiter, log_file)

        expected = "%s/expect_subset.fasta" % data_dir
        self.assertTrue(filecmp.cmp(out_fasta, expected, shallow=False))
        os.unlink(out_fasta)

        expected = "%s/expect_annotated_subset.csv" % data_dir
        self.assertTrue(filecmp.cmp(out_metadata, expected, shallow=False))
        os.unlink(out_metadata)