"""
This file will take the mapped_seq_id.csv file from id_reads
It will simplify its output into a new csv containing the following headers:
    ref_id, sequence_length, divergence, cluster_number, cluster_size

This file will be placed under results/id_reads/mapped_reads/simple_mapped_seq_id.csv
"""

import pandas as pd
from pprint import pprint


def calculate_cluster_sizes(data_frame):
    data_frame_clusters = data_frame["cluster"]
    cluster_counts = {}
    for cluster in data_frame_clusters:
        if cluster in cluster_counts:
            cluster_counts[cluster] += 1
        else:
            cluster_counts[cluster] = 1
    return cluster_counts


in_file = str(snakemake.input)
output = str(snakemake.output).split(" ")
divergence_threshold = snakemake.params.divergence_threshold
map_data_frame = pd.read_csv(in_file, header=0)

# create three pandas dataframes. One within bounds, one outside bounds, and one with NaN data
# store these in a list so we can iterate through each one and map it to its output file
divergence_data_frames = [
    map_data_frame.loc[map_data_frame["divergence"] <= divergence_threshold],
    map_data_frame.loc[map_data_frame["divergence"] > divergence_threshold],
    map_data_frame.loc[pd.isnull(map_data_frame["divergence"])]
]

for index, (data_frame, out_file) in enumerate(zip(divergence_data_frames, output)):
    print(f"Creating simple mapped sequence: {output}")
    cluster_size_dict = calculate_cluster_sizes(data_frame)
    ref_id = data_frame["ref_id"]
    sequence_length = data_frame["seq_length"]
    divergence = data_frame["divergence"]
    cluster_number = data_frame["cluster"]
    cluster_sizes = pd.Series([cluster_size_dict[item] for item in cluster_number])  # create a series matching the cluster size dictionary with cluster number pd.Series

    new_data = []
    for i, (ref_id, sequence_length, divergence, cluster_number, cluster_sizes) in enumerate(zip(ref_id, sequence_length, divergence, cluster_number, cluster_sizes)):
        new_data.append([ref_id, sequence_length, divergence, cluster_number, cluster_sizes])

    new_df = pd.DataFrame(new_data, columns=["ref_id", "seq_length", "divergence", "cluster", "cluster_size"])
    new_df.drop_duplicates(inplace=True)
    new_df.to_csv(out_file, index=False)

