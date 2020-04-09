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

def load_dataframe(metadata_file, filter_columns, where_columns):
    sep = ','
    if metadata_file.endswith('tsv'):
        sep = '\t'
    df = pd.read_csv(metadata_file, sep=sep)

    if where_columns:
        for pair in where_columns:
            column,regex = pair.split("=")
            df = find_column_with_regex(df, column, regex)

    df.rename(str.lower, axis='columns', inplace=True)

    if filter_columns:
        filter_columns = [s.lower() for s in filter_columns]
        df = df.loc[:, df.columns.isin(filter_columns)]
    return df

def add_data(new_dataframe, master_dataframe):
    column_intersection = [s for s in new_dataframe.columns if s in master_dataframe.columns]
    master_dataframe = master_dataframe.merge(new_dataframe, how='outer', on=column_intersection)
    return master_dataframe

def load_metadata(list_metadata_files, filter_columns, where_columns):
    master = None
    for metadata_file in list_metadata_files:
        if master is None:
            master = load_dataframe(metadata_file, filter_columns, where_columns)
        else:
            new_data = load_dataframe(metadata_file, filter_columns, where_columns)
            master = add_data(new_data, master)
    return master

def grouped_to_csv(grouped, index_columns, log_handle):
    index_columns.append("count")
    log_handle.write("%s\n" %','.join(index_columns))

    for name, group in grouped:
        if isinstance(name, str):
            v = [name]
        else:
            v = list(name)
        num_in_group = len(group.index)
        v.append(num_in_group)
        str_v = [str(s) for s in v]
        log_handle.write("%s\n" %','.join(str_v))

def get_groups(df, index_columns, log_handle=None):
    grouped = df.groupby(index_columns)
    if log_handle:
        grouped_to_csv(grouped, index_columns, log_handle)
    return grouped

def load_target_sample_sizes(target_file, index_columns):
    targets = {}
    df = pd.read_csv(target_file)
    df.set_index(index_columns, inplace=True)
    for name, row in df.iterrows():
        targets[name] = row["count"]
    return targets

def subsample_metadata(df, index_columns, sample_size, target_file, exclude_uk, log_handle=None):
    grouped = get_groups(df, index_columns, log_handle)
    targets = {}
    if target_file:
        targets = load_target_sample_sizes(target_file, index_columns)

    subsampled_indexes = []

    if exclude_uk:
        uk_mask = df['admin0'] == "UK"
        subsampled_indexes.extend(df[uk_mask].index)
        non_uk_df = df[~uk_mask]
        grouped = non_uk_df.groupby(index_columns)

    for name, group in grouped:
        num_in_group = len(group.index)
        target_size = sample_size
        if name in targets:
            target_size = targets[name]
        if num_in_group > target_size:
            subsampled_indexes.extend(group.sample(n=target_size, random_state=1).index)
        else:
            subsampled_indexes.extend(group.index)

    return df.iloc[subsampled_indexes]

def get_index_field_from_header(record, header_delimiter, index_field):
    if not index_field:
        return record.id

    fields = record.id.split(header_delim)

    if isinstance(index_field, int):
        assert index_field < len(fields)
        return fields[index_field]

    for field in fields:
        if field.startswith(index_field):
            return field[len(index_field)+1:]

    return record.id

