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
    def assertDataframeEqual(self, a, b, msg):
        try:
            pd_testing.assert_frame_equal(a, b)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def setUp(self):
        self.addTypeEqualityFunc(pd.DataFrame, self.assertDataframeEqual)

    def test_find_column_with_regex_column_exists_already(self):
        df = pd.DataFrame({'name': ['a','b','c','d'], "place": ['x','y','z','x1'],
                           "date": ['2020-04-01', '2020-04-05', '2020-03-29','2020-04-02']})
        column = 'name'
        regex =  '.lace'
        result = find_column_with_regex(df, column, regex)
        expect = df
        print(result)
        print(expect)
        self.assertEqual(result, expect)

    def test_find_column_with_regex_is_match(self):
        df = pd.DataFrame({'name': ['a','b','c','d'], "place": ['x','y','z','x1'],
                           "date": ['2020-04-01', '2020-04-05', '2020-03-29','2020-04-02']})
        column = 'Place'
        regex =  '.lace'
        result = find_column_with_regex(df, column, regex)
        expect = df
        expect['Place'] = expect['place']
        print(result)
        print(expect)
        self.assertEqual(result, expect)

    def test_find_column_with_regex_no_match(self):
        df = pd.DataFrame({'name': ['a','b','c','d'], "place": ['x','y','z','x1'],
                           "date": ['2020-04-01', '2020-04-05', '2020-03-29','2020-04-02']})
        column = 'Place'
        regex =  '.lant'
        result = find_column_with_regex(df, column, regex)
        expect = df
        print(result)
        print(expect)
        self.assertEqual(result, expect)

    def test_load_dataframe_csv(self):
        metadata_file = "%s/metadata.csv" %data_dir
        index_columns = None
        where_columns = None
        result = load_dataframe(metadata_file, index_columns, where_columns)
        expect = pd.DataFrame({'name': ['a','b','c','d'], "place": ['x','y','z','x1'],
                               "date": ['2020-04-01', '2020-04-05', '2020-03-29','2020-04-02']})
        print(result)
        print(expect)
        self.assertEqual(result, expect)

    def test_load_dataframe_tsv(self):
        metadata_file = "%s/metadata.tsv" %data_dir
        index_columns = None
        where_columns = None
        result = load_dataframe(metadata_file, index_columns, where_columns)
        expect = pd.DataFrame({'name': ['e','f','g','d'], "place": ['x','y','z','x1'],
                               "date": ['2020-04-01', '2020-04-05', '2020-03-29','2020-04-02'],
                               "blah": [1, 2, 3, 4]})
        print(result)
        print(expect)
        self.assertEqual(result, expect)

    def test_load_dataframe_tsv_where(self):
        metadata_file = "%s/metadata.tsv" %data_dir
        index_columns = None
        where_columns = ["Foo=[Bb]l.h"]
        result = load_dataframe(metadata_file, index_columns, where_columns)
        expect = pd.DataFrame({'name': ['e','f','g','d'], "place": ['x','y','z','x1'],
                               "date": ['2020-04-01', '2020-04-05', '2020-03-29','2020-04-02'],
                               "blah": [1, 2, 3, 4], "foo": [1, 2, 3, 4]})
        print(result)
        print(expect)
        self.assertEqual(result, expect)

    def test_load_dataframe_tsv_where_index(self):
        metadata_file = "%s/metadata.tsv" %data_dir
        index_columns = ["name", "date", "Foo"]
        where_columns = ["Foo=[Bb]l.h"]
        result = load_dataframe(metadata_file, index_columns, where_columns)
        expect = pd.DataFrame({'name': ['e','f','g','d'],
                               "date": ['2020-04-01', '2020-04-05', '2020-03-29','2020-04-02'],
                               "foo": [1, 2, 3, 4]})
        print(result)
        print(expect)
        self.assertEqual(result, expect)

    def test_load_metadata(self):
        list_metadata_files = ["%s/metadata.csv" %data_dir, "%s/metadata.tsv" %data_dir]
        index_columns = None
        where_columns = None
        result = load_metadata(list_metadata_files, index_columns, where_columns)
        expect = pd.DataFrame({'name': ['a','b','c','d','e','f','g'],
                               "place": ['x','y','z','x1','x','y','z'],
                               "date": ['2020-04-01', '2020-04-05', '2020-03-29','2020-04-02',
                                        '2020-04-01', '2020-04-05', '2020-03-29'],
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
        df = pd.DataFrame({'name': ['a', 'b', 'b', 'd'], "place": ['x', 'y', 'y', 'x1'],
                           "date": ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']})
        index_columns = ["name", "place"]
        target_file = None
        sample_size = 1
        select_by_max_column = None
        select_by_min_column = None
        exclude_uk = False
        log_handle = None
        result = subsample_metadata(df, index_columns, sample_size, target_file, select_by_max_column,
                                    select_by_min_column, exclude_uk, log_handle)
        result.reset_index(drop=True, inplace=True)
        expect1 = pd.DataFrame({'name': ['a', 'b', 'd'], "place": ['x', 'y', 'x1'],
                           "date": ['2020-04-01', '2020-03-29', '2020-04-02']})
        expect2 = pd.DataFrame({'name': ['a', 'b', 'd'], "place": ['x', 'y', 'x1'],
                                "date": ['2020-04-01', '2020-04-05', '2020-04-02']})
        if expect1.equals(result):
            self.assertEqual(result, expect1)
        else:
            self.assertEqual(result, expect2)

    def test_subsample_metadata_target_file(self):
        df = pd.DataFrame({'name': ['a', 'b', 'b', 'd'], "place": ['x', 'y', 'y', 'x1'],
                           "date": ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']})
        index_columns = ["name", "place"]
        target_file = "%s/target_counts.csv" %data_dir
        sample_size = 0
        select_by_max_column = None
        select_by_min_column = None
        exclude_uk = False
        log_handle = None
        result = subsample_metadata(df, index_columns, sample_size, target_file, select_by_max_column,
                                    select_by_min_column, exclude_uk, log_handle)
        result.reset_index(drop=True, inplace=True)
        expect1 = pd.DataFrame({'name': ['a', 'b', 'd'], "place": ['x', 'y', 'x1'],
                           "date": ['2020-04-01', '2020-03-29', '2020-04-02']})
        expect2 = pd.DataFrame({'name': ['a', 'b', 'd'], "place": ['x', 'y', 'x1'],
                                "date": ['2020-04-01', '2020-04-05', '2020-04-02']})
        if expect1.equals(result):
            self.assertEqual(result, expect1)
        else:
            self.assertEqual(result, expect2)

    def test_subsample_metadata_target_file_overrides_sample_size(self):
        df = pd.DataFrame({'name': ['a', 'b', 'b', 'd'], "place": ['x', 'y', 'y', 'x1'],
                           "date": ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']})
        index_columns = ["name", "place"]
        target_file = "%s/target_counts.csv" %data_dir
        sample_size = 4
        select_by_max_column = None
        select_by_min_column = None
        exclude_uk = False
        log_handle = None
        result = subsample_metadata(df, index_columns, sample_size, target_file, select_by_max_column,
                                    select_by_min_column, exclude_uk, log_handle)
        result.reset_index(drop=True, inplace=True)
        expect1 = pd.DataFrame({'name': ['a', 'b', 'd'], "place": ['x', 'y', 'x1'],
                           "date": ['2020-04-01', '2020-03-29', '2020-04-02']})
        expect2 = pd.DataFrame({'name': ['a', 'b', 'd'], "place": ['x', 'y', 'x1'],
                                "date": ['2020-04-01', '2020-04-05', '2020-04-02']})
        if expect1.equals(result):
            self.assertEqual(result, expect1)
        else:
            self.assertEqual(result, expect2)

    def test_get_index_field_from_header_no_index_field(self):
        record = SeqIO.read("%s/single.fasta" %data_dir, "fasta")
        index_field = None
        header_delimiter = '|'
        result = get_index_field_from_header(record, header_delimiter, index_field)
        expect = "COUNTRY/CODE/YEAR|ID||DATE|BLAH=XXX|CLOWN:YYY"
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

    def test_get_index_column_values_default(self):
        df = pd.DataFrame({'name': ['a', 'b', 'c', 'd'], "header": ['x', 'y', 'z', 'x1'],
                           "date": ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']})
        index_column = None
        df, result = get_index_column_values(df, index_column)
        expect = ['x', 'y', 'z', 'x1']
        self.assertEqual(result, expect)

    def test_get_index_column_values_default_no_header(self):
        df = pd.DataFrame({'name': ['a', 'b', 'c', 'd'], "place": ['x', 'y', 'z', 'x1'],
                           "date": ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']})
        index_column = None
        df, result = get_index_column_values(df, index_column)
        expect = ['a', 'b', 'c', 'd']
        self.assertEqual(result, expect)

    def test_get_index_column_values_int(self):
        df = pd.DataFrame({'name': ['a', 'b', 'c', 'd'], "place": ['x', 'y', 'z', 'x1'],
                           "date": ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']})
        index_column = [2]
        df, result = get_index_column_values(df, index_column)
        expect = ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']
        self.assertEqual(result, expect)

    def test_get_index_column_values_string(self):
        df = pd.DataFrame({'name': ['a', 'b', 'c', 'd'], "place": ['x', 'y', 'z', 'x1'],
                           "date": ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']})
        index_column = ["date"]
        df, result = get_index_column_values(df, index_column)
        expect = ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']
        self.assertEqual(result, expect)

    def test_get_index_column_values_int_too_big(self):
        df = pd.DataFrame({'name': ['a', 'b', 'c', 'd'], "place": ['x', 'y', 'z', 'x1'],
                           "date": ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']})
        index_column = [3]
        self.assertRaises(AssertionError, get_index_column_values, df, index_column)

    def test_get_index_column_values_string_not_column(self):
        df = pd.DataFrame({'name': ['a', 'b', 'c', 'd'], "place": ['x', 'y', 'z', 'x1'],
                           "date": ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']})
        index_column = ["id"]
        self.assertRaises(AssertionError, get_index_column_values, df, index_column)

    def test_get_index_column_values_column_has_duplicates(self):
        df = pd.DataFrame({'name': ['a', 'b', 'b', 'd'], "place": ['x', 'y', 'z', 'x1'],
                           "date": ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']})
        index_column = ["name"]
        #df, result = get_index_column_values(df, index_column)
        #self.assertEqual(len(df.index.values),2)
        self.assertRaises(AssertionError, get_index_column_values, df, index_column)

    def test_get_index_column_values_strings(self):
        df = pd.DataFrame({'name': ['a', 'b', 'c', 'd'], "place": ['x', 'y', 'z', 'x1'],
                           "date": ['2020-04-01', '2020-04-05', '2020-03-29', '2020-04-02']})
        index_column = ["name","date"]
        df, result = get_index_column_values(df, index_column)
        expect = ['a|2020-04-01', 'b|2020-04-05', 'c|2020-03-29', 'd|2020-04-02']
        self.assertEqual(result, expect)
