import pandas as pd
import os
from natsort import natsorted as natural_sort
from pprint import pprint


# split input/output files into lists
input_files = str(snakemake.input).split(" ")
output_files = str(snakemake.output).split(" ")

# iterate through input/output files, 'pairing' them together
for index, (in_file, out_file) in enumerate(zip(input_files, output_files)):
    print(f"Creating OTU table for: {os.path.basename(in_file)}")
    divergence_df = pd.read_csv(in_file, delimiter=",", header=0)

    """
    Use list comprehension for getting the following (example)
    bc1
    bc10
    bc9
    unclassified

    otu4
    otu9
    otu110

    we are using natural_sort() to sort numbers from smallest to largest (https://pypi.org/project/natsort/)

    Use a set to get unique barcodes/clusters
    """
    barcode_rows = natural_sort([f"bc{item}" for item in set(divergence_df["barcode"]) if item != "unclassified"])
    cluster_columns = natural_sort([f"otu{item}" for item in set(divergence_df["cluster"])])
    cluster_columns.insert(0, "")  # insert blank cell in header

    # Get count of duplicate barcodes/clusters in divergence CSV
    groupings = dict()
    for i, (barcode, cluster) in enumerate(zip(divergence_df["barcode"], divergence_df["cluster"])):
        # get barcode-otu mapping
        if "unclassified" not in str(barcode):
            mapping = f"bc{barcode}-otu{cluster}"
        else:
            continue  # continue to next item if barcode == unclassified

        # add mapping to grouping
        if mapping in groupings:
            # add one to existing group
            groupings[mapping] += 1
        else:
            # make new group
            groupings[mapping] = 1

    # we are trying to create following table (example)
    """
            otu1    otu2    otu3    .   .   .
    bc1     1       4       2       .   .   .
    bc2     9       8       4       .   .   .
    bc3     4       1       6       .   .   .
    bc4     6       7       3       .   .   .
    .       .       .       .       .   .   .
    .       .       .       .       .   .   .
    .       .       .       .       .   .   .        
    """

    # create a matrix to load our values into
    otu_data = []
    for col in range(len(barcode_rows)):
        temp_list = []
        for row in range(len(cluster_columns) - 1):
            temp_list.append(0)
        otu_data.append(temp_list)

    # prepend each otu row with its barcode
    for i, (barcode_row, otu_row) in enumerate(zip(barcode_rows, otu_data)):
        otu_row.insert(0, barcode_row)

    # iterate through groupings
    for i, key in enumerate(groupings):
        barcode_index = 0

        # get the current group barcode/cluster we are working on
        group_barcode = key.split("-")[0]
        group_cluster = key.split("-")[1]

        # iterate through barcodes in this divergence table
        for current_barcode in barcode_rows:

            # set cluster_index to 1 to skip the 'barcode' header (as this will never match an otu value)
            cluster_index = 1

            # only iterate through clusters if barcodes match (should speed up matching a bit)
            if current_barcode == group_barcode:

                # iterate through clusters
                for current_cluster in cluster_columns[1:]:

                    # if grouping barcode/cluster match the barcode/cluster index, update index
                    if (current_barcode == group_barcode) and (current_cluster == group_cluster):
                        otu_data[barcode_index][cluster_index] = groupings[key]

                    cluster_index += 1  # increment cluster/barcode to next index
            barcode_index += 1

    # create data frame
    otu_data_frame = pd.DataFrame(data=otu_data, columns=cluster_columns)
    otu_data_frame.to_csv(path_or_buf=out_file, index=False)
