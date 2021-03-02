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

def find_column_with_regex(df, column, regex):
    regex = re.compile(regex)
    for original_column in df.columns:
        match = re.search(regex, original_column)
        if match and column in df.columns:
            print("Update column", column, "with non-na values from column", original_column)
            df[column].update(df[original_column])
            return df
        elif match:
            print("New column", column, "with non-na values from column", original_column)
            df[column] = df[original_column]
    return df

def find_field_with_regex(header, regex):
    regex = re.compile(regex)
    match = re.search(regex, header)
    if match:
        return match.group()
    else:
        return ""

def load_dataframe(metadata_file, filter_columns, where_columns, omit_columns=False):
    sep = ','
    na_values = ["None"]
    if metadata_file.endswith('tsv'):
        sep = '\t'
    try:
        df = pd.read_csv(metadata_file, sep=sep, na_values=na_values)
    except:
        df = pd.read_csv(metadata_file, sep=sep, na_values=na_values, encoding='utf-8')

    if where_columns:
        for pair in where_columns:
            column,regex = pair.split("=")
            df = find_column_with_regex(df, column, regex)

    if 'unnamed: 0' in df.columns:
        df.drop(columns=['unnamed: 0'], inplace=True)

    if filter_columns:
        df = filter_by_omit_columns(df)
        df = df.loc[:, df.columns.isin(filter_columns)]
        for column in [c for c in filter_columns if c not in df.columns.values]:
            df[column] = None
    else:
        df.dropna(how='all', axis='columns', inplace=True)

    return df

def add_data(new_dataframe, master_dataframe):
    column_intersection = [s for s in new_dataframe.columns if s in master_dataframe.columns]
    master_dataframe = master_dataframe.merge(new_dataframe, how='outer', on=column_intersection)
    return master_dataframe

def get_header_id(record, where_field):
    if where_field is None:
        return record.id
    else:
        field, regex = where_field.split("=")
        return find_field_with_regex(record.description, regex)

def filter_by_omit_columns(df):
    drop_indexes = []
    for column in df.columns.values:
        if "omit" in column.lower():
            drop_indexes.extend(df.index[df[column] == True].tolist())
    df = df.drop(drop_indexes)
    return df

def load_metadata(list_metadata_files, filter_columns=None, where_columns=None, index_column=None):
    master = None
    for metadata_file in list_metadata_files:
        if master is None:
            master = Metadata(metadata_file, where_columns, filter_columns, index_column)
        else:
            new_data = Metadata(metadata_file, where_columns, filter_columns, index_column)
            master.add_data(new_data)
    return master

def load_metadata_df(list_metadata_files, filter_columns=None, where_columns=None, index_column=None):
    master = None
    for metadata_file in list_metadata_files:
        if master is None:
            master = load_dataframe(metadata_file, filter_columns, where_columns)
        else:
            new_data = load_dataframe(metadata_file, filter_columns, where_columns)
            master = add_data(new_data, master)
    if filter_columns:
        master = master[filter_columns]
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

def add_subsample_omit_column(df, non_omitted_df, subsampled_df):
    if "subsample_omit" not in df.columns:
        df["subsample_omit"] = False
    non_omitted_df_index_values = non_omitted_df.index.values
    subsampled_df_index_values = subsampled_df.index.values
    for i in non_omitted_df_index_values:
        if i not in subsampled_df_index_values:
            df.loc[i,"subsample_omit"] = True
    return

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

def get_index_column_values(df, index_columns, header_delimiter='|'):
    print(index_columns, type(index_columns))
    if isinstance(index_columns, str):
    	index_columns = [index_columns]
    if not index_columns or len(index_columns) == 0 or index_columns == "":
        if "sequence_name" in df.columns:
            index_columns = ["sequence_name"]
        else:
            index_columns = [0]

    str_index_columns = []
    for column in index_columns:
        if isinstance(column, int):
            assert column < len(df.columns.values)
            column = df.columns[column]
        print(column, df.columns.values)
        assert column in list(df.columns.values)
        str_index_columns.append(column)

    column_values = []
    for i,row in df.iterrows():
        column_values.append(header_delimiter.join([str(row[c]) for c in str_index_columns]))
    index_column_name = "|".join(str_index_columns)
    df.loc[:,index_column_name] = column_values
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
