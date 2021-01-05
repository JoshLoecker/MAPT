# -*- coding: utf-8 -*-
"""
"""

# define inputs
import pandas as pd
import numpy as np
import os
import pysam  # for loading .sam files
import sys
import re  # string stuff

# import seaborn as sns
# from matplotlib import pyplot as plt


def get_tag(tag_list, tag):
    """
    Extract the value of a SAM tag

    input:
        tag_list: str pd.Series containing tags (usually dicts)
        tag: str identifier of the tag to extract

    returns:
        the value of the tag, probably as a string but maybe in correct type.
    """
    has_tag = tag in str(tag_list)

    if has_tag:
        out = [i.split(':')[2] for i in tag_list if tag in str(i)]
        return out[0]


def sam_to_DataFrame(samfile, primary_only=True):
    """
    Load a file using pysam, and convert it to a pd.DataFrame object.

    input:
        samfile: str path to samfile to load
        primary_only: logical only keep primary alignments?

    returns:
        pd.DataFrame
    """
    out = [i.to_dict() for i in samfile.fetch(until_eof=True)]

    out = pd.DataFrame(out)

    num_cols = ['flag',
                'ref_pos',
                'map_quality',
                'next_ref_pos',
                'length',
                'qual']
    out[num_cols] = out[num_cols].apply(pd.to_numeric, errors='coerce')

    if primary_only:
        tt = out['tags'].apply(get_tag, tag='tp') == 'P'
        out = out[tt]

    return out


def fc(DataFrame, n=0):
    """
    View pd.DataFrame.iloc[n, ]
    """
    out = DataFrame.iloc[n, ]
    return out


def sub_ls(pattern, replacement, string_list):
    """
    iteratively call re.sub

    inputs:
        pattern: str to search for
        replacement: str
        string_list: list (of strings)

    returns:
        list (of strings)
    """
    out = [re.sub(pattern, replacement, i) for i in string_list]
    return out


def collapse_tags(x, sep='; '):
    """
    Concatenate all items in x into a single string.

    input:
        x: accepts any input, but most useful for something that's iterable.
        sep: separation between items in x

    returns:
        str
    """

    def _is_iterable(obj):
        try:
            iter(obj)
        except Exception:
            rr = False
        else:
            rr = True
        return rr

    if not _is_iterable(x):
        out = x
    else:
        out = sep.join(x)
    return out


def extract_cluster(x, style='isONclust'):
    """
    extracts cluster number from x

    input:
        x: pd.Series (strings) containing sampe info
        style: str. Specify the program that created x. spoa or isONclust.

    returns:
        cluster number from x.
    """
    if style == 'isONclust':
        out = x.split(':')[1]
        out = out.split('_')[0]
    elif style == 'spoa':
        out = x.split('_')[1]
    else:
        raise ValueError("Style must be 'isONclust' (default) or 'spoa'")

    return int(out)


# load reads-by-clusters info
clusters = str(snakemake.input.clustering) + str("/final_clusters.tsv")
clusters = pd.read_csv(clusters, sep='\t', engine="python")
clusters.columns = ('cluster', 'read_id')

# format
find_barcode = 'barcode='
clusters['barcode'] = [i.split(find_barcode)[1] for i in clusters['read_id']]
clusters['barcode'] = [i.split('_')[0] for i in clusters['barcode']]
clusters['read_id'] = [i.split('_')[0] for i in clusters['read_id']]

# load minimap output (aligned consensus sequences)
mapping = snakemake.input.minimap
mapping = pysam.AlignmentFile(mapping, 'r')
mapping = sam_to_DataFrame(mapping, primary_only=True)

# reformat alignment output
mapping['cluster'] = mapping['name'].apply(extract_cluster, style='spoa')

has_de = mapping['tags'].str.contains('de')

mapping['divergence'] = mapping['tags'].apply(get_tag, tag='de')  # defined above
mapping['divergence_type'] = 'gap-compressed'
mapping['score'] = mapping['tags'].apply(get_tag, tag='AS')
mapping['mismatches'] = mapping['tags'].apply(get_tag, tag='NM')
mapping['chimeric'] = mapping['tags'].apply(get_tag, tag='SA')
mapping['chimeric'] = mapping['chimeric']
mapping['seq_length'] = mapping['seq'].apply(len)

mapping['ref_id'] = [i.split('_1')[0] for i in mapping['ref_name']]

mapping.to_csv(snakemake.output.minimap_output)

# merge mapping and clusters info
out = pd.merge(clusters, mapping, how='left', on='cluster')

replace_names = {
    'read_id'    : 'seq_id',
    'flag'       : 'bit_flag',
    'map_quality': 'quality'}

for i in replace_names.keys():
    out.columns = sub_ls(i, replace_names[i], out.columns)

out['aligner'] = 'ison_minimap'

# subset columns
keep = ['seq_id',
        'ref_id',
        'bit_flag',
        'seq_length',
        'quality',
        'divergence',
        'divergence_type',
        'score',
        'mismatches',
        'chimeric',
        'aligner',
        'barcode',
        'cluster',
        'tags']
out = out[keep]

out['barcode'] = sub_ls('barcode', '', out['barcode'])

out['tags'] = out['tags'].apply(collapse_tags)

# export alignment to csv for each sequence table for downstream.
# This is the file we're really interested in.
file_out = snakemake.output.mapped_seq_id_csv
out.to_csv(file_out, index=False)

# summarize the minimap output and export to csv
keep = ['cluster',
        'ref_id',
        'bit_flag',
        'divergence',
        'seq_length',
        'tags']
out_uniques = out.groupby(by=keep).size().reset_index(name='n_reads')

file_out = snakemake.output.mapped_consensus_csv
out_uniques.to_csv(file_out, index=False)
