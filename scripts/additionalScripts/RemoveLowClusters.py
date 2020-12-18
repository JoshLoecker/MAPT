import os
import shutil
print("Starting. . .")
input = "/Users/joshl/Downloads/temp_results/isONclust/cluster_fastq"
output = "/Users/joshl/Downloads/temp_results/temp/"
process_done = "/Users/joshl/Downloads/temp_results/temp/RemoveLowClustersDone"

read_range = 18
reads = [0 for x in range(0, read_range)]

for file in sorted(os.listdir(input)):
    path = os.path.join(input, file)
    file_lines = open(path, 'r').readlines()
    total_lines = len(file_lines)
    total_reads = total_lines / 4

    if total_reads <= 2:
        shutil.move(src=path, dst=output)


# 'touch' output file
open(output, "w").close()
