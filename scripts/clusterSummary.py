import pandas as pd

def calculate_cluster_sizes(clusters):
    """
    This function will create a dictionary of clusters
    Key: Cluster number
    Value: Frequency of cluster within data

    :param clusters: This is a Pandas series of the "cluster" header
    :return: a dictionary of cluster numbers as keys and its frequency as values
    """
    counts = {}
    for clust in clusters:
        if clust in counts:
            counts[clust] += 1
        else:
            counts[clust] = 1
    return counts


output_files = str(snakemake.output).split(" ")
map_data_frame = pd.read_csv(str(snakemake.input), header=0)

divergence_data_frames = [
    map_data_frame.loc[map_data_frame["divergence"] <= divergence_threshold],
    map_data_frame.loc[map_data_frame["divergence"] > divergence_threshold],
    map_data_frame.loc[pd.isnull(map_data_frame["divergence"])]
]

for index, (data_frame, out_file) in enumerate(zip(divergence_data_frames, output_files)):
    cluster_numbers = data_frame["cluster"]
    cluster_counts_dict = calculate_cluster_sizes(cluster_numbers)
    cluster_sizes = pd.Series([cluster_counts_dict[item] for item in cluster_numbers])  # create a series matching the cluster size dictionary with cluster number pd.Series

    new_data = []
    for i, (number, size) in enumerate(zip(cluster_numbers, cluster_sizes)):
        new_data.append([number, size])

    new_df = pd.DataFrame(new_data, columns=["cluster", "size"])
    new_df.to_csv(out_file, index=False)

