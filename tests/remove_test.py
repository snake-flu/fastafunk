import os
import unittest
import filecmp

from fastafunk.remove import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'remove')

class TestRemove(unittest.TestCase):
    def test_run_remove(self):
        in_fasta = ["%s/test_remove.fasta" %data_dir]
        in_metadata = ["%s/test_filter.tsv" %data_dir]
        out_fasta = "%s/tmp.filtered.fasta" %data_dir
        log_file = "%s/tmp.filtered.fasta.log" %data_dir
        expected = "%s/filtered.fasta" %data_dir
        remove_fasta(in_fasta, in_metadata, out_fasta, log_file)
        self.assertTrue(filecmp.cmp(out_fasta, expected, shallow=False))
        os.unlink(out_fasta)

        self.assertTrue(filecmp.cmp(out_fasta + ".log", expected + ".log", shallow=False))
        os.unlink(out_fasta + ".log")