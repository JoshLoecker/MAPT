from pprint import pprint
import re
input_files = str(snakemake.input).split(" ")



for file in input_files:
    if re.search("barcode..\.merged.fastq", file):
        print("CLASSIFIED")
    elif re.search("unclassified.merged.fastq", file):
        print("UN CLASSIFIED")
