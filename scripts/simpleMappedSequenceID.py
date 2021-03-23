"""
This file will take the mapped_seq_id.csv file from id_reads
It will simplify its output into a new csv containing the following headers:
    ref_id, sequence_length, divergence, cluster_number, cluster_size

This file will be placed under results/id_reads/mapped_reads/simple_mapped_seq_id.csv
"""

import pandas as pd
from pprint import pprint

map_data_frame = pd.read_csv(str(snakemake.input), header=0)
output_files = str(snakemake.output).split(" ")
divergence_threshold = snakemake.params.divergence_threshold

# create three pandas dataframes. One within bounds, one outside bounds, and one with NaN data
# store these in a list so we can iterate through each one and map it to its output file
divergence_data_frames = [
    map_data_frame.loc[map_data_frame["divergence"] <= divergence_threshold],
    map_data_frame.loc[map_data_frame["divergence"] > divergence_threshold],
    map_data_frame.loc[pd.isnull(map_data_frame["divergence"])]
]

for index, (data_frame, out_file) in enumerate(zip(divergence_data_frames, output_files)):
    print(f"Creating simple mapped sequence: {out_file}")
    ref_id = data_frame["ref_id"]
    sequence_length = data_frame["seq_length"]
    divergence = data_frame["divergence"]
    cluster_number = data_frame["cluster"]
    otu = [f"otu{i}" for i in cluster_number] #PME

    new_data = []
    for i, (ref_id, sequence_length, divergence, cluster_number, otu) in enumerate(zip(ref_id, sequence_length, divergence, cluster_number, otu)):
        new_data.append([ref_id, sequence_length, divergence, cluster_number, otu])

    new_df = pd.DataFrame(new_data, columns=["ref_id", "seq_length", "divergence", "cluster", "otu"])
    new_df.drop_duplicates(inplace=True)
    new_df.to_csv(out_file, index=False)

