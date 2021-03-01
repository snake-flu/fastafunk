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

from fastafunk.utils import *

this_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
data_dir = os.path.join(this_dir, 'tests', 'data', 'utils')

class TestUtils(unittest.TestCase):
    def assertDataFrameEqual(self, a, b, msg):
        try:
            pd_testing.assert_frame_equal(a, b)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def setUp(self):
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataFrameEqual)

    def test_find_column_with_regex_column_exists_already(self):
        metadata_file = "%s/metadata.csv" %data_dir
        metadata = Metadata(metadata_file)
        column = 'name'
        regex =  '.lace'
        metadata.find_column_with_regex(column, regex)
        result = "%s/tmp.find_column_with_regex_column_exists_already.csv" %data_dir
        result_handle = open(result,'w')
        metadata.to_csv(result_handle)
        result_handle.close()
        expect = "%s/metadata_overwritten_name.csv" %data_dir
        self.assertTrue(filecmp.cmp(result, expect, shallow=False))
        os.unlink(result)

    def test_find_column_with_regex_is_match(self):
        metadata_file = "%s/metadata.csv" %data_dir
        metadata = Metadata(metadata_file)
        column = 'Place'
        regex =  '.lace'
        metadata.find_column_with_regex(column, regex)
        result = "%s/tmp.find_column_with_regex_is_match.csv" %data_dir
        result_handle = open(result,'w')
        metadata.to_csv(result_handle)
        result_handle.close()
        expect = "%s/metadata_with_Place.csv" %data_dir
        self.assertTrue(filecmp.cmp(result, expect, shallow=False))
        os.unlink(result)

    def test_find_column_with_regex_no_match(self):
        metadata_file = "%s/metadata.csv" %data_dir
        metadata = Metadata(metadata_file)
        column = 'Place'
        regex =  '.lant'
        metadata.find_column_with_regex(column, regex)
        result = "%s/tmp.find_column_with_regex_no_match.csv" %data_dir
        result_handle = open(result,'w')
        metadata.to_csv(result_handle)
        result_handle.close()
        expect = "%s/metadata.csv" %data_dir
        self.assertTrue(filecmp.cmp(result, expect, shallow=False))
        os.unlink(result)

    def test_find_field_with_regex_no_match(self):
        header = ">hCoV-19/place/code/2020||country|date"
        regex = "EPI_[\w_]+"
        result = find_field_with_regex(header, regex)
        expect = ""
        self.assertEqual(result, expect)

    def test_find_field_with_regex_match(self):
        header = ">EPI_ISL_422243|hCoV-19/place/code/2020||country|date"
        regex = "EPI_[\w_]+"
        result = find_field_with_regex(header, regex)
        expect = "EPI_ISL_422243"
        self.assertEqual(result, expect)

    def test_find_field_with_regex_no_match2(self):
        header = ">EPI_ISL_422243|place/code/2020||country|date"
        regex = "hCo[vV]-19/\w+/[\w_-]+/\w+"
        result = find_field_with_regex(header, regex)
        expect = ""
        self.assertEqual(result, expect)

    def test_find_field_with_regex_match2(self):
        header = ">EPI_ISL_422243|hCoV-19/place/code/2020||country|date"
        regex = "hCo[vV]-19/\w+/[\w_-]+/\w+"
        result = find_field_with_regex(header, regex)
        expect = "hCoV-19/place/code/2020"
        self.assertEqual(result, expect)

    def test_find_field_with_regex_match3(self):
        header = ">hCoV-19/place/code/2020||country|date"
        regex = "hCo[vV]-19/\w+/[\w_-]+/\w+"
        result = find_field_with_regex(header, regex)
        expect = "hCoV-19/place/code/2020"
        self.assertEqual(result, expect)

    def test_load_metadata_csv(self):
        metadata_file = "%s/metadata.csv" %data_dir
        index_columns = None
        where_columns = None
        result = load_metadata([metadata_file], index_columns, where_columns)
        result = pd.DataFrame(result.rows)
        expect = pd.DataFrame({'name': ['a','b','c','d'], "place": ['x','y','z','x1'],
                               "date": ['2020-04-01', '2020-04-05', '2020-03-29','2020-04-02']})
        print(result)
        print(expect)
        self.assertEqual(result, expect)

    def test_load_metadata_tsv(self):
        metadata_file = "%s/metadata.tsv" %data_dir
        index_columns = None
        where_columns = None
        result = load_metadata([metadata_file], index_columns, where_columns)
        result = pd.DataFrame(result.rows)
        expect = pd.DataFrame({'Name': ['e','f','g','d'], "Place": ['x','y','z','x1'],
                               "Date": ['2020-04-01', '2020-04-05', '2020-03-29','2020-04-02'],
                               "Blah": ['1', '2', '3', '4']})
        print(result)
        print(expect)
        self.assertEqual(result, expect)

    def test_load_metadata_tsv_where(self):
        metadata_file = "%s/metadata.tsv" %data_dir
        index_columns = None
        where_columns = ["Foo=[Bb]l.h"]
        result = load_metadata([metadata_file], index_columns, where_columns)
        result = pd.DataFrame(result.rows)
        expect = pd.DataFrame({'Name': ['e','f','g','d'], "Place": ['x','y','z','x1'],
                               "Date": ['2020-04-01', '2020-04-05', '2020-03-29','2020-04-02'],
                               "Blah": ['1', '2', '3', '4'], "Foo": ['1', '2', '3', '4']})
        print(result)
        print(expect)
        self.assertEqual(result, expect)

    def test_load_metadata_tsv_where_index(self):
        metadata_file = "%s/metadata.tsv" %data_dir
        index_columns = ["Name", "Date", "Foo"]
        where_columns = ["Foo=[Bb]l.h"]
        result = load_metadata([metadata_file], index_columns, where_columns)
        result = pd.DataFrame(result.rows)
        expect = pd.DataFrame({'Name': ['e','f','g','d'],
                               "Date": ['2020-04-01', '2020-04-05', '2020-03-29','2020-04-02'],
                               "Foo": ['1', '2', '3', '4']})
        print(result)
        print(expect)
        self.assertEqual(result, expect)

    def test_filter_by_omit_columns(self):
        metadata_file = "%s/metadata_with_omit.csv" %data_dir
        metadata = Metadata(metadata_file)
        metadata.filter_by_omit_columns()
        result = "%s/tmp.filter_by_omit_columns.csv" %data_dir
        result_handle = open(result,'w')
        metadata.to_csv(result_handle)
        result_handle.close()
        expect = "%s/metadata_filtered_by_omit.csv" %data_dir
        self.assertTrue(filecmp.cmp(result, expect, shallow=False))
        os.unlink(result)

    def test_load_metadata(self):
        list_metadata_files = ["%s/metadata.csv" %data_dir, "%s/metadata.tsv" %data_dir]
        index_columns = None
        where_columns = ["name=Name","place=Place", "date=Date", "blah=Blah"]
        result = load_metadata(list_metadata_files, index_columns, where_columns)
        result = pd.DataFrame(result.rows)
        result.Blah = result.Blah.astype(float)
        result.blah = result.blah.astype(float)
        expect = pd.DataFrame({'name': ['a','b','c','d','e','f','g'],
                               "place": ['x','y','z','x1','x','y','z'],
                               "date": ['2020-04-01', '2020-04-05', '2020-03-29','2020-04-02',
                                        '2020-04-01', '2020-04-05', '2020-03-29'],
                               'Name': [None, None, None,'d','e','f','g'],
                               "Place": [None, None, None,'x1','x','y','z'],
                               "Date": [None, None, None,'2020-04-02',
                                        '2020-04-01', '2020-04-05', '2020-03-29'],
                               "Blah": [None, None, None, 4, 1, 2, 3],
                               "blah": [None, None, None, 4, 1, 2, 3]})
        print(result)
        print(expect)
        self.assertEqual(result, expect)

    def test_grouped_to_csv(self):
        df = pd.DataFrame({'name': ['a', 'b', 'b', 'd'], "place": ['x', 'y', 'y', 'x1'],
                           "date": ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']})
        index_columns = ["name","place"]
        grouped = df.groupby(index_columns)
        result = "%s/tmp.group_counts.csv" %data_dir
        result_handle = open(result,'w')
        grouped_to_csv(grouped, index_columns, result_handle)
        result_handle.close()
        expect = "%s/group_counts.csv" %data_dir
        self.assertTrue(filecmp.cmp(result, expect, shallow=False))
        os.unlink(result)

    def test_load_target_sample_sizes(self):
        target_file = "%s/group_counts.csv" %data_dir
        index_columns = ["name","place"]
        targets = load_target_sample_sizes(target_file, index_columns)
        expect = {('a', 'x'): 1, ('b', 'y'): 2, ('d', 'x1'): 1}
        self.assertEqual(targets, expect)

    def test_get_subsample_indexes_max(self):
        group = df = pd.DataFrame({'name': ['a', 'b', 'b', 'd'], "place": ['x', 'y', 'y', 'x1'],
                                   "missing": [3,5,7,3], "length": [10,2,11,12],
                                   "date": ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']})
        target_size = 1
        select_by_max_column = "length"
        select_by_min_column = None
        result = get_subsample_indexes(group, target_size, select_by_max_column, select_by_min_column)
        expect = [3]
        self.assertEqual(result, expect)

    def test_get_subsample_indexes_min(self):
        group = df = pd.DataFrame({'name': ['a', 'b', 'b', 'd'], "place": ['x', 'y', 'y', 'x1'],
                                   "missing": [3,5,7,3], "length": [10,2,10,10],
                                   "date": ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']})
        target_size = 1
        select_by_max_column = None
        select_by_min_column = "missing"
        result = get_subsample_indexes(group, target_size, select_by_max_column, select_by_min_column)
        expect1 = [0]
        expect2 = [3]
        if result == expect1:
            self.assertEqual(result, expect1)
        else:
            self.assertEqual(result, expect2)

    def test_get_subsample_indexes_na(self):
        group = df = pd.DataFrame({'name': ['a', 'b', 'b', 'd'], "place": ['x', 'y', 'y', 'x1'],
                                   "missing": [3,5,None,3], "length": [10,2,None,10],
                                   "date": ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']})
        target_size = 3
        select_by_max_column = None
        select_by_min_column = "missing"
        result = get_subsample_indexes(group, target_size, select_by_max_column, select_by_min_column)
        expect = [0,1,3]
        self.assertEqual(result, expect)

    def test_get_subsample_indexes_na_max_target_size(self):
        group = df = pd.DataFrame({'name': ['a', 'b', 'b', 'd'], "place": ['x', 'y', 'y', 'x1'],
                                   "missing": [3,5,None,3], "length": [10,2,None,10],
                                   "date": ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']})
        target_size = 4
        select_by_max_column = None
        select_by_min_column = "missing"
        result = get_subsample_indexes(group, target_size, select_by_max_column, select_by_min_column)
        expect = [0,1,3]
        self.assertEqual(result, expect)

    def test_get_subsample_indexes_max_target_size(self):
        group = df = pd.DataFrame({'name': ['a', 'b', 'b', 'd'], "place": ['x', 'y', 'y', 'x1'],
                                   "missing": [3,5,None,3], "length": [10,2,None,10],
                                   "date": ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']})
        target_size = 4
        select_by_max_column = None
        select_by_min_column = None
        result = get_subsample_indexes(group, target_size, select_by_max_column, select_by_min_column)
        expect = [0,1,2,3]
        self.assertEqual(result, expect)

    def test_subsample_metadata_sample_size(self):
        index_column = "name"
        target_file = None
        sample_size = 1
        select_by_max_column = None
        select_by_min_column = None
        exclude_uk = False
        log_handle = None

        metadata_file = "%s/metadata_to_subsample.csv" %data_dir
        metadata = load_metadata_df([metadata_file])
        subsampled_metadata = subsample_metadata(metadata, index_column, sample_size, target_file, select_by_max_column,
                           select_by_min_column, exclude_uk, log_handle)
        add_subsample_omit_column(metadata, metadata, subsampled_metadata)
        result = "%s/tmp.subsample_metadata_sample_size.csv" %data_dir
        result_handle = open(result,'w')
        metadata.to_csv(result_handle, index=False)
        result_handle.close()
        expect1 = "%s/metadata_subsampled1.csv" %data_dir
        expect2 = "%s/metadata_subsampled2.csv" %data_dir
        if filecmp.cmp(result, expect1, shallow=False):
            self.assertTrue(filecmp.cmp(result, expect1, shallow=False))
        else:
            self.assertTrue(filecmp.cmp(result, expect2, shallow=False))
        os.unlink(result)

    def test_subsample_metadata_target_file(self):
        index_column = "place"
        target_file = "%s/target_counts.csv" %data_dir
        sample_size = 0
        select_by_max_column = None
        select_by_min_column = None
        exclude_uk = False
        log_handle = None

        metadata_file = "%s/metadata_to_subsample.csv" %data_dir
        metadata = load_metadata_df([metadata_file])
        subsampled_metadata = subsample_metadata(metadata, index_column, sample_size, target_file, select_by_max_column,
                                                 select_by_min_column, exclude_uk, log_handle)
        add_subsample_omit_column(metadata, metadata, subsampled_metadata)
        result = "%s/tmp.subsample_metadata_target_file.csv" %data_dir
        result_handle = open(result,'w')
        metadata.to_csv(result_handle, index=False)
        result_handle.close()
        expect1 = "%s/metadata_subsampled1.csv" %data_dir
        expect2 = "%s/metadata_subsampled2.csv" %data_dir
        if filecmp.cmp(result, expect1, shallow=False):
            self.assertTrue(filecmp.cmp(result, expect1, shallow=False))
        else:
            self.assertTrue(filecmp.cmp(result, expect2, shallow=False))
        os.unlink(result)

    def test_subsample_metadata_target_file_overrides_sample_size(self):
        index_column = "place"
        target_file = "%s/target_counts.csv" %data_dir
        sample_size = 4
        select_by_max_column = None
        select_by_min_column = None
        exclude_uk = False
        log_handle = None

        metadata_file = "%s/metadata_to_subsample.csv" %data_dir
        metadata = load_metadata_df([metadata_file])
        subsampled_metadata = subsample_metadata(metadata, index_column, sample_size, target_file, select_by_max_column,
                                                 select_by_min_column, exclude_uk, log_handle)
        add_subsample_omit_column(metadata, metadata, subsampled_metadata)
        result = "%s/tmp.subsample_metadata_target_file_overrides_sample_size.csv" %data_dir
        result_handle = open(result,'w')
        metadata.to_csv(result_handle, index=False)
        result_handle.close()
        expect1 = "%s/metadata_subsampled1.csv" %data_dir
        expect2 = "%s/metadata_subsampled2.csv" %data_dir
        if filecmp.cmp(result, expect1, shallow=False):
            self.assertTrue(filecmp.cmp(result, expect1, shallow=False))
        else:
            self.assertTrue(filecmp.cmp(result, expect2, shallow=False))
        os.unlink(result)

    def test_get_index_field_from_header_no_index_field(self):
        record = SeqIO.read("%s/single.fasta" %data_dir, "fasta")
        index_field = None
        header_delimiter = '|'
        result = get_index_field_from_header(record, header_delimiter, index_field)
        expect = "COUNTRY/CODE/YEAR|ID||DATE|BLAH=XXX|CLOWN:YYY foo=bar parrot:clown"
        self.assertEqual(result, expect)

    def test_get_index_field_from_header_int_index_field(self):
        record = SeqIO.read("%s/single.fasta" %data_dir, "fasta")
        index_field = 0
        header_delimiter = '|'
        result = get_index_field_from_header(record, header_delimiter, index_field)
        expect = "COUNTRY/CODE/YEAR"
        self.assertEqual(result, expect)

    def test_get_index_field_from_header_int_index_field1(self):
        record = SeqIO.read("%s/single.fasta" %data_dir, "fasta")
        index_field = 1
        header_delimiter = '|'
        result = get_index_field_from_header(record, header_delimiter, index_field)
        expect = "ID"
        self.assertEqual(result, expect)

    def test_get_index_field_from_header_int_index_field_empty(self):
        record = SeqIO.read("%s/single.fasta" %data_dir, "fasta")
        index_field = 2
        header_delimiter = '|'
        result = get_index_field_from_header(record, header_delimiter, index_field)
        expect = ""
        self.assertEqual(result, expect)

    def test_get_index_field_from_header_str_index_field(self):
        record = SeqIO.read("%s/single.fasta" %data_dir, "fasta")
        index_field = "BLAH"
        header_delimiter = '|'
        result = get_index_field_from_header(record, header_delimiter, index_field)
        expect = "XXX"
        self.assertEqual(result, expect)

    def test_get_index_field_from_header_str_index_field1(self):
        record = SeqIO.read("%s/single.fasta" %data_dir, "fasta")
        index_field = "CLOWN"
        header_delimiter = '|'
        result = get_index_field_from_header(record, header_delimiter, index_field)
        expect = "YYY"
        self.assertEqual(result, expect)

    def test_get_index_field_from_header_str_index_field_header_delim(self):
        record = SeqIO.read("%s/single.fasta" %data_dir, "fasta")
        index_field = "foo"
        header_delimiter = ' '
        result = get_index_field_from_header(record, header_delimiter, index_field)
        expect = "bar"
        self.assertEqual(result, expect)

    def test_get_index_field_from_header_str_index_field_header_delim1(self):
        record = SeqIO.read("%s/single.fasta" %data_dir, "fasta")
        index_field = "parrot"
        header_delimiter = ' '
        result = get_index_field_from_header(record, header_delimiter, index_field)
        expect = "clown"
        self.assertEqual(result, expect)

    def test_get_index_column_values_default_sequence_name(self):
        metadata_file = "%s/metadata_with_sequence_name.csv" %data_dir
        metadata = Metadata(metadata_file)
        result = metadata.get_index_column_values()
        expect = ['x', 'y', 'z', 'x1']
        self.assertEqual(result, expect)

    def test_get_index_column_values_default_0(self):
        metadata_file = "%s/metadata.csv" %data_dir
        metadata = Metadata(metadata_file)
        result = metadata.get_index_column_values()
        expect = ['a', 'b', 'c', 'd']
        self.assertEqual(result, expect)

    def test_get_index_column_values_int(self):
        metadata_file = "%s/metadata.csv" %data_dir
        metadata = Metadata(metadata_file,index=2)
        result = metadata.get_index_column_values()
        expect = ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']
        self.assertEqual(result, expect)

    def test_get_index_column_values_string(self):
        metadata_file = "%s/metadata.csv" %data_dir
        metadata = Metadata(metadata_file,index="date")
        result = metadata.get_index_column_values()
        expect = ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']
        self.assertEqual(result, expect)

    def test_set_index_column_int_too_big(self):
        metadata_file = "%s/metadata.csv" %data_dir
        self.assertRaises(IndexError, Metadata, metadata_file,index=3)

    def test_set_index_column_string_not_column(self):
        metadata_file = "%s/metadata.csv" %data_dir
        self.assertRaises(SystemExit, Metadata, metadata_file,index="id")

