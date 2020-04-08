import sys
import csv
import pandas as pd
import re

def get_log_handle(log_file, out_fasta):
    if log_file:
        log_handle = open(log_file,"w")
    elif out_fasta:
        log_handle = sys.stdout
    else:
        log_handle = sys.stderr
    return log_handle

def get_out_handle(out_fasta):
    if not out_fasta:
        out_handle = sys.stdout
    else:
        out_handle = open(out_fasta,"w")
    return out_handle

def get_in_handle(in_fasta):
    if not in_fasta:
        in_handle = sys.stdin
    else:
        in_handle = open(in_fasta,"r")
    return in_handle

def close_handle(handle):
    if handle:
        handle.close()

def metadata_to_dict(list_metadata_files):
    metadata_dictionary = {}

    for metadata_file in list_metadata_files:
        sep = ','
        if metadata_file.endswith('tsv'):
            sep = '\t'

        with open(metadata_file) as csv_handle:
            csv_reader = csv.reader(csv_handle, delimiter=sep)
            for row in csv_reader:
                metadata_dictionary[row[0]] = row[1:]

    return metadata_dictionary

def find_column_with_regex(df, column, regex):
    if column in df.columns:
        return df
    regex = re.compile(regex)
    for original_column in df.columns:
        match = re.search(regex, original_column)
        print(original_column, regex, match)
        if match:
            df[column] = df[original_column]
            return df
    return df

def load_dataframe(metadata_file, index_columns, where_columns):
    sep = ','
    if metadata_file.endswith('tsv'):
        sep = '\t'
    df = pd.read_csv(metadata_file, sep=sep)

    print(df)
    if where_columns:
        for pair in where_columns:
            column,regex = pair.split("=")
            print(pair, column, regex)
            df = find_column_with_regex(df, column, regex)
            print(df)

    df.rename(str.lower, axis='columns', inplace=True)
    print(df)

    if index_columns:
        index_columns = [s.lower() for s in index_columns]
        df = df.loc[:, df.columns.isin(index_columns)]
    print(df)
    return df

def load_metadata(list_metadata_files, index_columns, where_columns):
    master = None
    master_key = None

    for metadata_file in list_metadata_files:
        if master is None:
            master = load_dataframe(metadata_file, index_columns, where_columns)
        else:
            new_data = load_dataframe(metadata_file, index_columns, where_columns)
            column_intersection = [s for s in new_data.columns if s in master.columns]
            master = master.merge(new_data, how='outer', on=column_intersection)
    return master

