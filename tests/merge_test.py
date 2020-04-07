import os
import unittest
import filecmp
import glob

from fastafunk.merge import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'merge')

class TestMerge(unittest.TestCase):
    def test_run_merge_with_log(self):
        in_fasta = glob.glob("%s/*.f*a" %data_dir)
        in_metadata = "%s/metadata.csv" %data_dir
        out_fasta = "%s/tmp.merged_fasta" %data_dir
        log_file = "%s/tmp.merged_fasta.log" %data_dir
        expected = "%s/merged_fasta" %data_dir
        merge_fasta(in_fasta, in_metadata, out_fasta, log_file)
        self.assertTrue(filecmp.cmp(out_fasta, expected, shallow=False))
        os.unlink(out_fasta)
        os.unlink(log_file)

    def test_run_merge_no_log(self):
        in_fasta = glob.glob("%s/*.f*a" % data_dir)
        in_metadata = "%s/metadata.csv" % data_dir
        out_fasta = "%s/tmp.merged_fasta" % data_dir
        log_file = None
        expected = "%s/merged_fasta" % data_dir
        merge_fasta(in_fasta, in_metadata, out_fasta, log_file)
        self.assertTrue(filecmp.cmp(out_fasta, expected, shallow=False))
        os.unlink(out_fasta)

    def test_run_merge_no_log_no_out_fasta(self):
        in_fasta = glob.glob("%s/*.f*a" % data_dir)
        in_metadata = "%s/metadata.csv" % data_dir
        out_fasta = ""
        log_file = None
        expected = "%s/merged_fasta" % data_dir
        merge_fasta(in_fasta, in_metadata, out_fasta, log_file)
        self.assertTrue(filecmp.cmp(out_fasta, expected, shallow=False))
        os.unlink(out_fasta)
