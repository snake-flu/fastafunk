import csv
import re
import sys
from datetime import datetime

def clean_dict(d, column_names=None):
    to_delete = []
    for key in d.keys():
        if key == '' or "unnamed" in key:
            to_delete.append(key)
        elif column_names is not None and key not in column_names:
            to_delete.append(key)
    for key in to_delete:
        del d[key]
    return d

class MetadataReader:

    def __init__(self, metadata_file, where_columns=None, filter_columns=None, index=None):
        self.columns = []
        self.where_column_dict = {}
        self.rows = []
        self.omit_rows = []
        self.index = None
        self.reader = None
        self.sep = ","
        self.file = metadata_file
        self.handle = open(metadata_file, "r")

        self.load_from_file(metadata_file, where_columns, filter_columns, index)


    def get_index(self, index):
        if index is None:
            if "sequence_name" in self.columns:
                self.index = "sequence_name"
            else:
                self.index = self.columns[0]
        elif index in self.columns:
            self.index = index
        elif isinstance(index, int):
            self.index = self.columns[index]
        else:
            sys.exit("Error, index %s is not a column" %index)

    def get_columns(self,where_columns,filter_columns):
        self.columns = [c for c in self.reader.fieldnames if c != '' and "unnamed" not in c]
        if where_columns:
            for pair in where_columns:
                column,regex = pair.split("=")
                if column not in self.where_column_dict:
                    self.where_column_dict[column] = []
                regex = re.compile(regex)
                for original_column in self.columns:
                    match = re.search(regex, original_column)
                    if match:
                        self.where_column_dict[column].append(original_column)
                if column not in self.columns:
                    self.columns.append(column)
        if filter_columns:
            self.columns = filter_columns

    def get_rows(self):
        omit_columns = [ c for c in self.reader.fieldnames if "omit" in c.lower() or c in ["duplicate", "why_excluded"]]
        for row in self.reader:
            omit = False
            for column in omit_columns:
                if row[column] not in ["False", False, None, "None", ""]:
                    omit = True
                    continue
            if omit:
                self.omit_rows.append(row[self.index])
            else:
                self.rows.append(row[self.index])

    def get_reader(self):
        self.reader = csv.DictReader(self.handle, delimiter=self.sep)

    def load_from_file(self, metadata_file, where_columns=None, filter_columns=None, index=None):
        if metadata_file.endswith('tsv'):
            self.sep = '\t'
        self.get_reader()
        self.get_columns(where_columns, filter_columns)
        self.get_index(index)
        self.get_rows()
        self.close()
        self.handle = open(self.file, "r")
        self.get_reader()

    def add_columns(self, new_columns):
        self.columns.extend([c for c in new_columns if c not in self.columns])

    def clean_row(self, row_dict, new_data_dict=None):
        clean_values = {}
        for column in self.columns:
            value = None
            if column in self.where_column_dict:
                for where_column in self.where_column_dict[column]:
                    if row_dict[where_column] not in [None, "None", ""]:
                        value = row_dict[where_column]
            elif column in row_dict:
                value = row_dict[column]
            if new_data_dict is not None and column in new_data_dict:
                value = new_data_dict[column]
            clean_values[column] = value
        return clean_values

    def restrict(self, sequence_list):
        self.omit_rows.extend(self.rows)
        self.rows = sequence_list

    def to_csv(self, out_handle, header=True, include_omitted=False, new_data_dict=None):
        self.get_reader()
        writer = csv.DictWriter(out_handle, fieldnames = self.columns, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        if header:
            writer.writeheader()
        for row in self.reader:
            if include_omitted or row[self.index] not in self.omit_rows:
                if new_data_dict is not None and row[self.index] in new_data_dict:
                    writer.writerow(self.clean_row(row, new_data_dict[row[self.index]]))
                else:
                    writer.writerow(self.clean_row(row, None))

    def close(self):
        self.handle.close()





