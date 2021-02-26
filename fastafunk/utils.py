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
import dendropy
from fastafunk.metadata import *


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

def trees_to_taxa(list_tree_files):
    """
    extract all taxa labels in newick or nexus format trees
    labels is a list of taxon labels
    """
    labels = []
    for treefile in list_tree_files:
        try:
            tree = dendropy.Tree.get(path = treefile,
                                     schema = 'newick',
                                     preserve_underscores=True)
        except:
            try:
                tree = dendropy.Tree.get(path = treefile,
                                         schema = 'nexus',
                                         preserve_underscores=True)
            except:
                sys.exit('tree does not seem to be in Newick or Nexus format')

        labels = labels + tree.taxon_namespace.labels()

    return(labels)

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

def load_metadata(list_metadata_files, filter_columns=None, where_columns=None, index_column=None):
    master = None
    for metadata_file in list_metadata_files:
        if master is None:
            master = Metadata(metadata_file, where_columns, filter_columns, index_column)
        else:
            new_data = Metadata(metadata_file, where_columns, filter_columns, index_column)
            master.add_data(new_data)
    return master

def load_new_metadata(list_metadata_files, date_column, filter_columns=None, where_columns=None, index_column=None):
    master = None
    date = None

    for metadata_file in list_metadata_files:
        if master is None:
            master = Metadata(metadata_file, where_columns, filter_columns, index_column)
            date = master.get_newest_date(date_column)
            if date is None:
                master.fill_date_where_missing(date_column)
                date = master.get_newest_date(date_column)
        else:
            new_data = Metadata(metadata_file, where_columns, filter_columns, index_column)
            new_date = new_data.get_newest_date(date_column)
            if new_date is None:
                master.add_data(new_data)
                master.fill_date_where_missing(date_column)
                master.subset_by_min(date_column,date)
            elif new_date > date:
                master.subset_by_min(date_column,date)
                date = new_date
            elif new_date < date:
                mmaster.subset_by_min(date_column,new_date)
    if len(master.rows) == 0:
        sys.exit("Check date column %s exists in at least one file")
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

def subsample_metadata(metadata, group_columns, sample_size, target_file, select_by_max_column, select_by_min_column,
                       exclude_uk, log_handle=None):
    non_omit_rows = metadata.get_omit_rows(inverse=True)
    df = pd.DataFrame(metadata.rows)
    if select_by_max_column:
        df.select_by_max_column.dropna(inplace=True)
        data_types_dict = {select_by_max_column: int}
        df = df.astype(data_types_dict)
    if select_by_min_column:
        df.select_by_min_column.dropna(inplace=True)
        data_types_dict = {select_by_min_column: int}
        df = df.astype(data_types_dict)
    df.iloc[non_omit_rows, :]
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

    metadata.add_subsample_column(subsampled_indexes)
    subsampled = df.loc[subsampled_indexes]
    return subsampled[metadata.index].values

def get_index_field_from_header(record, header_delimiter, index_field):
    if index_field is None or index_field == "" or index_field == []:
        return record.description

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

def get_cov_id(record):
    id_strings = record.id.split('|')[0].split('/')
    if len(id_strings) > 2 and id_strings[0].startswith("hCo"):
        id_string = id_strings[2]
    elif len(id_strings) > 1:
        id_string = id_strings[1]
    else:
        id_string = ""
    return id_string

def restrict_dataframe(df, column_name, values_to_keep):
    mask = [not x[0] for x in df[column_name].isin(values_to_keep).values.tolist()]
    drop_indexes = df.index[mask].tolist()
    df.drop(drop_indexes, inplace=True)
    return df

def clean_dict(d, column_names=None):
    to_delete = []
    for key in d.keys():
        if key == '':
            to_delete.append(key)
        elif "unnamed" in key:
            to_delete.append(key)
        elif column_names is not None and key not in column_names:
            to_delete.append(key)
    for key in to_delete:
        del d[key]
    return d
