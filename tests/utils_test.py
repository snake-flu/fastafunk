import os
import unittest
import filecmp
import pandas as pd
import pandas.testing as pd_testing

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
                               "blah": ['NaN','NaN','NaN', 4.0, 1.0, 2.0, 3.0]})
        print(result)
        print(expect)
        self.assertEqual(result, expect)