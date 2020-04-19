"""
Name: utils.py
Author: Rachel Colquhoun
Date: 08 April 2020
Description: Utility functions used within Fastafunk

This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import sys
import csv
import pandas as pd
import numpy as np
import re
import datetime

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
        in_handle = open(in_fasta)
    return in_handle

def close_handle(handle):
    if handle:
        handle.close()

def fix_header_string(header):
    return header\
        .replace("hCoV-19/","")\
        .replace("hCov-19/","")\
        .replace(" ","_")

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
        if match:
            df[column] = df[original_column]
            return df
    return df

def find_field_with_regex(header, regex):
    regex = re.compile(regex)
    match = re.search(regex, header)
    if match:
        return match.group()
    else:
        return ""

def get_header_id(record, where_field):
    if where_field is None:
        return record.id
    else:
        field, regex = where_field.split("=")
        return find_field_with_regex(record.description, regex)

def load_dataframe(metadata_file, filter_columns, where_columns):
    sep = ','
    if metadata_file.endswith('tsv'):
        sep = '\t'
    try:
        df = pd.read_csv(metadata_file, sep=sep)
    except:
        df = pd.read_csv(metadata_file, sep=sep, encoding='utf-8')

    if where_columns:
        for pair in where_columns:
            column,regex = pair.split("=")
            df = find_column_with_regex(df, column, regex)

    df.rename(str.lower, axis='columns', inplace=True)

    if filter_columns:
        filter_columns = [s.lower() for s in filter_columns]
        df = df.loc[:, df.columns.isin(filter_columns)]

    if 'unnamed: 0' in df.columns:
        df.drop(columns=['unnamed: 0'], inplace=True)
    df.dropna(how='all', axis='columns', inplace=True)
    return df

def add_data(new_dataframe, master_dataframe):
    column_intersection = [s for s in new_dataframe.columns if s in master_dataframe.columns]
    master_dataframe = master_dataframe.merge(new_dataframe, how='outer', on=column_intersection)
    return master_dataframe

def add_empty_columns(new_columns, master_dataframe):
    column_difference = [s for s in new_columns if s not in master_dataframe.columns]
    for column in column_difference:
        master_dataframe[column] = ""
    return master_dataframe

def filter_by_omit_columns(df):
    drop_indexes = []
    for column in df.columns.values:
        if "omit" in column.lower():
            drop_indexes.extend(df.index[df[column] == True].tolist())
    df = df.drop(drop_indexes)
    return df

def load_metadata(list_metadata_files, filter_columns, where_columns):
    master = None
    for metadata_file in list_metadata_files:
        if master is None:
            master = load_dataframe(metadata_file, filter_columns, where_columns)
        else:
            new_data = load_dataframe(metadata_file, filter_columns, where_columns)
            master = add_data(new_data, master)
    return master

def get_newest_date(df, column):
    if column not in df.columns.values:
        return None
    df[column] = pd.to_datetime(df[column])
    newest_date = df[column].max()
    return newest_date

def fill_date_where_missing(df,column):
    date = datetime.date.today()
    if column not in df.columns.values:
        df[column] = date
    df[column].fillna(date, inplace=True)
    return df

def load_new_metadata(list_metadata_files, date_column, filter_columns=None, where_columns=None):
    master = None
    date = None

    for metadata_file in list_metadata_files:
        if master is None:
            master = load_dataframe(metadata_file, filter_columns, where_columns)
            date = get_newest_date(master, date_column)
        else:
            new_data = load_dataframe(metadata_file, filter_columns, where_columns)
            new_date = get_newest_date(new_data, date_column)
            if new_date is None:
                master = add_data(new_data, master)
                master = fill_date_where_missing(master, date_column)
            elif date is None:
                master = add_data(master, new_data)
                master = fill_date_where_missing(master, date_column)
            elif new_date > date:
                master = new_data.loc[pd.to_datetime(new_data[date_column]) > date]
                date = new_date
            elif new_date < date:
                master = master.loc[pd.to_datetime(master[date_column]) > new_date]

    return master

def grouped_to_csv(grouped, group_columns, log_handle):
    group_columns.append("count")
    log_handle.write("%s\n" %','.join(group_columns))

    for name, group in grouped:
        if isinstance(name, str):
            v = [name]
        else:
            v = list(name)
        num_in_group = len(group.index)
        v.append(num_in_group)
        str_v = [str(s) for s in v]
        log_handle.write("%s\n" %','.join(str_v))

def get_groups(df, group_columns, log_handle=None):
    grouped = df.groupby(group_columns)
    if log_handle:
        grouped_to_csv(grouped, group_columns, log_handle)
    return grouped

def load_target_sample_sizes(target_file, group_columns):
    targets = {}
    df = pd.read_csv(target_file)
    df.set_index(group_columns, inplace=True)
    for name, row in df.iterrows():
        targets[name] = row["count"]
    return targets

def get_subsample_indexes(group, target_size, select_by_max_column, select_by_min_column):
    subsampled_indexes = []
    idx = 0

    if select_by_max_column and select_by_max_column in group.columns:
        while len(subsampled_indexes) < target_size and not np.isnan(idx):
            idx = group[select_by_max_column].idxmax(axis=0, skipna=True)
            if not np.isnan(idx):
                subsampled_indexes.append(idx)
                group.drop(idx, inplace=True)
        return sorted(subsampled_indexes)

    elif select_by_min_column and select_by_min_column in group.columns:
        while len(subsampled_indexes) < target_size and not np.isnan(idx):
            idx = group[select_by_min_column].idxmin(axis=0, skipna=True)
            if not np.isnan(idx):
                subsampled_indexes.append(idx)
                group.drop(idx, inplace=True)
        return sorted(subsampled_indexes)

    else:
        return sorted(list(group.sample(n=target_size, random_state=1).index))

def subsample_metadata(df, group_columns, sample_size, target_file, select_by_max_column, select_by_min_column,
                       exclude_uk, log_handle=None):
    grouped = get_groups(df, group_columns, log_handle)
    targets = {}
    if target_file:
        targets = load_target_sample_sizes(target_file, group_columns)

    subsampled_indexes = []

    if exclude_uk:
        uk_mask = df['admin0'] == "UK"
        subsampled_indexes.extend(df[uk_mask].index)
        non_uk_df = df[~uk_mask]
        grouped = non_uk_df.groupby(group_columns)

    for name, group in grouped:
        num_in_group = len(group.index)
        target_size = sample_size
        if name in targets:
            target_size = targets[name]
        if num_in_group > target_size:
            subsampled_indexes.extend(get_subsample_indexes(group, target_size, select_by_max_column,
                                                            select_by_min_column))
        else:
            subsampled_indexes.extend(group.index)

    subsampled_indexes.sort()
    return df.loc[subsampled_indexes]

def add_subsample_omit_column(df, subsampled_df):
    if "subsample_omit" not in df.columns:
        df["subsample_omit"] = True
    for i in subsampled_df.index.values:
        df.loc[i,"subsample_omit"] = False
    return

def get_index_field_from_header(record, header_delimiter, index_field):
    if index_field is None or index_field == "":
        return record.id

    if header_delimiter == " ":
        fields = record.description.split()
    else:
        fields = record.id.split(header_delimiter)


    if isinstance(index_field, int):
        assert index_field < len(fields)
        return fields[index_field]

    for field in fields:
        if field.startswith(index_field):
            return field[len(index_field)+1:]

    return record.id

def get_index_column_values(df, index_columns, header_delimiter='|'):
    if not index_columns or len(index_columns) == 0:
        if "sequence_name" in df.columns:
            index_columns = ["sequence_name"]
        else:
            index_columns = [0]

    str_index_columns = []
    for column in index_columns:
        if isinstance(column, int):
            assert column < len(df.columns.values)
            column = df.columns[column]
        assert column in df.columns.values
        str_index_columns.append(column)

    column_values = []
    for i,row in df.iterrows():
        column_values.append(header_delimiter.join([str(row[c]) for c in str_index_columns]))
    df.loc[:,"sequence_name"] = column_values
    #bad_headers = df[df["sequence_name"].duplicated()]["sequence_name"].index.values
    #df.drop(bad_headers, inplace=True)
    #if df["sequence_name"].duplicated().any():
    #    print(df.loc[df["sequence_name"].duplicated(),"sequence_name"])
    #assert not df["sequence_name"].duplicated().any()

    return df, column_values

def get_cov_id(record):
    id_strings = record.id.split('|')[0].split('/')
    if len(id_strings) > 2 and id_strings[0].startswith("hCo"):
        id_string = id_strings[2]
    elif len(id_strings) > 1:
        id_string = id_strings[1]
    else:
        id_string = ""
    return id_string
