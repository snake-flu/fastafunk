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

class Metadata:

    def __init__(self, metadata_file=None, where_columns=None, filter_columns=None, index=None, metadata_dict=None):
        self.columns = []
        self.rows = []
        self.index = None
        if metadata_file:
            self.load_from_file(metadata_file, where_columns, filter_columns, index)
        elif metadata_dict:
            self.load_from_dict(metadata_dict, index)

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

    def load_from_file(self,metadata_file, where_columns=None, filter_columns=None, index=None):
        with open(metadata_file,"r") as f:
            sep = ','
            if metadata_file.endswith('tsv'):
                sep = '\t'
            reader = csv.DictReader(f, delimiter=sep)
            self.columns = [c for c in reader.fieldnames if c != '' and "unnamed" not in c]
            self.rows = [clean_dict(r) for r in reader]

        self.apply_where_columns(where_columns)
        self.apply_filter_columns(filter_columns)
        self.get_index(index)

    def load_from_dict(self, data_dict, index=None):
        self.columns = list(data_dict.keys())
        self.get_index(index)

        for i in range(len(data_dict[index])):
            d = {}
            for c in data_dict:
                d[c] = data_dict[c][i]
            self.rows.append(d)

    def add_column(self, column, template=None):
        new_column = False
        if column not in self.columns:
            self.columns.append(column)
            new_column = True
        for row in self.rows:
            if template and row[template] not in [None, "None", ""]:
                row[column] = row[template]
            elif new_column:
                row[column] = None

    def add_null_columns(self, columns):
        self.columns.extend(columns)
        for row in self.rows:
            for column in columns:
                row[column] = None

    def find_column_with_regex(self, column, regex):
        regex = re.compile(regex)
        for original_column in self.columns:
            match = re.search(regex, original_column)
            if match:
                self.add_column(column, original_column)

    def apply_where_columns(self, where_columns=None):
        if where_columns:
            for pair in where_columns:
                column,regex = pair.split("=")
                self.find_column_with_regex(column, regex)

    def remove_columns(self, columns_to_remove):
        for key in columns_to_remove:
            for row in self.rows:
                del row[key]
            self.columns.remove(key)

    def apply_filter_columns(self, filter_columns=None):
        if filter_columns:
            self.remove_columns([c for c in self.columns if c not in filter_columns])

    def index_rows(self):
        indexed_rows = {}
        count = 0
        for row in self.rows:
            indexed_rows[row[self.index]] = count
            count += 1
        return indexed_rows

    def add_data(self, new_data):
        new_columns = [x for x in new_data.columns if x not in self.columns]
        self.add_null_columns(new_columns)
        #old_columns = [x for x in self.columns if x not in new_data.columns]
        #new_data.add_null_columns(old_columns)

        indexed_rows = self.index_rows()
        for row in new_data.rows:
            if row[new_data.index] in indexed_rows:
                self.rows[indexed_rows[row[new_data.index]]].update(row)
            else:
                self.rows.append(row)

    def get_newest_date(self, column):
        newest_date = None
        if column not in self.columns:
            return None
        for row in self.rows:
            row_date = datetime.fromisoformat(row[column])
            if newest_date is None or row_date > newest_date:
                newest_date = row_date
        return newest_date

    def fill_date_where_missing(column):
        date = datetime.date.today()
        if column not in self.columns:
            self.columns.append(column)
        for row in self.rows:
            if row[column] in [None, "None", ""]:
                row[column] = date

    def subset_by_min(self,column, min_value):
        if column not in self.columns:
            sys.exit("Error, column %s is not a column" %column)
        del_rows = []
        for i,row in enumerate(self.rows):
            if row[column] <= min_value:
                del_rows.append(i)
        for i in reversed(del_rows):
            self.rows.pop(i)

    def get_omit_rows(self, inverse=False):
        omit_columns = [column for column in self.columns if "omit" in column.lower()]
        del_rows = []
        for i,row in enumerate(self.rows):
            omit = False
            for column in omit_columns:
                if row[column] == "True":
                    omit = True
                    break
            if not inverse and omit:
                del_rows.append(i)
            elif inverse and not omit:
                del_rows.append(i)
        return del_rows

    def filter_by_omit_columns(self):
        del_rows = self.get_omit_rows()
        for i in reversed(del_rows):
            self.rows.pop(i)

    def add_subsample_column(self, keep_rows):
        self.add_column("subsample_omit")
        for i in reversed(range(len(self.rows))):
            if i not in keep_rows:
                self.rows[i]["subsample_omit"] = True

    def restrict(self, column, values_to_keep):
        if column not in self.columns:
            sys.exit("Error, column %s is not a column" %column)
        del_rows = []
        for i,row in enumerate(self.rows):
            if row[column] not in values_to_keep:
                del_rows.append(i)
        for i in reversed(del_rows):
            self.rows.pop(i)

    def get_index_column_values(self):
        values = []
        for row in self.rows:
            values.append(row[self.index])
        return values

    def to_csv(self, out_handle):
        writer = csv.DictWriter(out_handle, fieldnames = self.columns, delimiter=",", quotechar='\"', quoting=csv.QUOTE_MINIMAL, dialect = "unix")
        writer.writeheader()
        for row in self.rows:
            writer.writerow(row)





