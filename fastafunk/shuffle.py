"""
Name: shuffle.py
Author: Rachel Colquhoun
Date: 19 March 2021
Description: This module shuffles rows of a metadata file.

This file is part of Fastafunk (https://github.com/cov-ert/fastafunk).
Copyright 2020 Xiaoyu Yu (xiaoyu.yu@ed.ac.uk) & Rachel Colquhoun (rachel.colquhoun@ed.ac.uk).
"""

import random


def shuffle(in_metadata, out_metadata):
    with open(in_metadata, 'r') as in_csv, open(out_metadata, 'w') as out_csv:
        first_line = in_csv.readline()
        lines = in_csv.readlines()[1:]
        random.shuffle(lines)
        out_csv.write(first_line)
        out_csv.writelines(lines)
