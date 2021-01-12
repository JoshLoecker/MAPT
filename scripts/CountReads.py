"""
This file was created because we are counting reads at multiple spots in the pipeline
As a result, it is easier to make one function more robust than to copy and paste the same code multiple times
"""

# import libraries
import csv
import os
from collections import OrderedDict
import gzip

# iterate through and open each file in the input stream
total_reads_dict = {}
for file_path in snakemake.input:

    # get the total number of lines in the file
    # if the file is a .fastq.gz file, we will get a `UnicodeDecodeError`
    # use gzip instead in this case
    try:
        file_lines = open(file_path, 'r').readlines()
    except UnicodeDecodeError:
        file_lines = gzip.open(file_path, "rt").readlines()
    total_lines = len(file_lines)

    # we know that fastq files have four lines per read, and fasta files have two lines per read
    # we will divide the file by the appropriate number
    if ".fastq" in file_path:
        number_of_reads = int(total_lines / 4)
    elif "alignment.summary.csv" in file_path:
        number_of_reads = total_lines
    elif ".minimap.sam" in file_path:
        number_of_reads = 0
        # iterate through each line in the file
        for line in file_lines:
            # only count reads that do not start with "@"
            if line[0] != "@":
                number_of_reads += 1
    else:
        number_of_reads = 0

    # collect the barcode number from the file path
    file_name = os.path.basename(file_path)  # get the file name
    barcode_number = file_name.split(".")[0]  # get the barcode number

    # now add the barcode number with the total number of reads to our dictionary
    total_reads_dict[barcode_number] = number_of_reads

# sort the reads so it is easy to find them in the resulting file
# https://docs.python.org/2/library/collections.html#collections.OrderedDict
total_reads_dict = OrderedDict( sorted(total_reads_dict.items(), key=lambda x: x[0]) )

# now write the reads to our csv file
with open(snakemake.output[0], 'w') as output_stream:
    csv_writer = csv.writer(output_stream, delimiter="\t")

    # write a simple header row
    csv_writer.writerow(["barcode", "reads", "process"])

    # write our barcode number and reads in barcode in the format of: barcode## READS
    for barcode in total_reads_dict.keys():
        csv_writer.writerow([barcode, total_reads_dict[barcode], snakemake.params.process])
