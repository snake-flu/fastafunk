import os
import unittest

from fastafunk.stats import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'utils')

class TestStats(unittest.TestCase):
    def test_remove_terminal_gaps_dashs(self):
        sequence = "-----AGC-T---"
        result = remove_terminal_gaps(sequence)
        expect = "AGC-T"
        self.assertEqual(expect,result)

    def test_remove_terminal_gaps_ns(self):
        sequence = "NNNNNNAGC-T---"
        result = remove_terminal_gaps(sequence)
        expect = "AGC-T"
        self.assertEqual(expect,result)

    def test_get_length_unaligned_exclude_terminal_gaps(self):
        sequence = "NNNNNNAGC-T---"
        result = get_length(sequence)
        expect = 4
        self.assertEqual(expect, result)

    def test_get_length_unaligned_include_terminal_gaps(self):
        sequence = "NNNNNNAGC-T---"
        result = get_length(sequence, exclude_terminal_gaps=False)
        expect = 10
        self.assertEqual(expect, result)

    def test_get_length_aligned_exclude_terminal_gaps(self):
        sequence = "NNNNNNAGC-T---"
        result = get_length(sequence, unaligned=False)
        expect = 5
        self.assertEqual(expect, result)

    def test_get_length_aligned_include_terminal_gaps(self):
        sequence = "NNNNNNAGC-T---"
        result = get_length(sequence, unaligned=False, exclude_terminal_gaps=False)
        expect = 14
        self.assertEqual(expect, result)

    def test_get_number_missing_bases_unaligned_exclude_terminal_gaps(self):
        sequence = "NNNNNNAGC-T---"
        result = get_number_missing_bases(sequence)
        expect = 0
        self.assertEqual(expect, result)

    def test_get_number_missing_bases_unaligned_include_terminal_gaps(self):
        sequence = "NNNNNNAGC-T---"
        result = get_number_missing_bases(sequence, exclude_terminal_gaps=False)
        expect = 6
        self.assertEqual(expect, result)

    def test_get_number_missing_bases_aligned_exclude_terminal_gaps(self):
        sequence = "NNNNNNAGCN-T---"
        result = get_number_missing_bases(sequence, unaligned=False)
        expect = 1
        self.assertEqual(expect, result)

    def test_get_number_missing_bases_aligned_include_terminal_gaps(self):
        sequence = "NNNNNNAGCN-T---"
        result = get_number_missing_bases(sequence, unaligned=False, exclude_terminal_gaps=False)
        expect = 7
        self.assertEqual(expect, result)

    def test_get_proportion_gaps_unaligned_exclude_terminal_gaps(self):
        sequence = "NNNNNNAGCN-T---"
        result = get_proportion_gaps(sequence)
        expect = 1.0/5
        self.assertEqual(expect, result)

    def test_get_proportion_gaps_unaligned_include_terminal_gaps(self):
        sequence = "NNNNNNAGCN-T---"
        result = get_proportion_gaps(sequence, exclude_terminal_gaps=False)
        expect = 7.0/11
        self.assertEqual(expect, result)

    def test_get_proportion_gaps_aligned_exclude_terminal_gaps(self):
        sequence = "NNNNNNAGCN-T---"
        result = get_proportion_gaps(sequence, unaligned=False)
        expect = 1.0/6
        self.assertEqual(expect, result)

    def test_get_proportion_gaps_aligned_include_terminal_gaps(self):
        sequence = "NNNNNNAGCN-T---"
        result = get_proportion_gaps(sequence, unaligned=False, exclude_terminal_gaps=False)
        expect = 7.0/15
        self.assertEqual(expect, result)

    def test_get_stat(self):
        sequence = "NNNNNNAGCN-T---"
        result = get_stat("gaps", sequence, unaligned=False, exclude_terminal_gaps=False)
        expect = 7.0/15
        self.assertEqual(expect, result)