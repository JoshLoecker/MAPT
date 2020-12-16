import os
import subprocess
from pathlib import Path
import glob
from multiprocessing import cpu_count
from pprint import pprint

configfile: "config.yaml"

# these environment variables are required by this pipeline
envvars:
    "alignment_path",
    "num_basecallers",
    "threads_per_caller",
    "basecall_configuration",
    "barcode_kit",
    "cutadapt_trim_error_rate",
    "cutadapt_trim_three_prime_adapter",
    "cutadapt_trim_five_prime_adapter",
    "nanofilt_filtering_min",
    "nanofilt_filtering_max"


def return_barcode_numbers(path: str):
    """
    This function will return a list of barcode numbers under the directory passed in
    Args:
        path: The directory that should be searched
    Returns: A list of strings containing directory names
    """
    barcode_numbers = set()
    for item in os.scandir(path):
        item = item.name
        if "barcode" in item:
            barcode_numbers.add(item)
        elif "unclassified" in item:
            barcode_numbers.add(item)
    return barcode_numbers
def barcode_merge_files(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    barcodes = set()  # a set is like a list but only stores unique values
    for folder in os.listdir(checkpoint_output):
        full_path = os.path.join(checkpoint_output, folder)
        if Path(full_path).is_dir():
            barcodes.add(folder)

    return_merged_files = [config['results'] + "barcode/" + barcode + ".merged.fastq" for barcode in barcodes]
    return return_merged_files
def nanoqc_basecall_data(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config['results'] + "visuals/nanoqc/basecall/"
def nanoqc_barcode_classified(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config['results'] + "visuals/nanoqc/barcode/classified"
def nanoqc_barcode_unclassified(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config['results'] + "visuals/nanoqc/barcode/unclassified"
def cutadapt(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return expand(
        config['results'] + "cutadapt/{barcode}.cutadapt.fastq",
        barcode=return_barcode_numbers(checkpoint_output))
def filtering(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return expand(config['results'] + "filter/{barcode}.filter.fastq",
                  barcode=return_barcode_numbers(checkpoint_output))
def isONclust_pipeline(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config['results'] + "isONclust/pipeline/"
def isONclust_cluster_fastq(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config['results'] + "isONclust/cluster_fastq/"
def IsoCon(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config['results'] + "isoCon/"
def guppy_aligner(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return expand(
        config['results'] + "alignment/guppy/sam_files/{barcode}.guppy.sam",
               barcode=return_barcode_numbers(checkpoint_output))
def minimap_aligner_from_filtering(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return expand(
        config['results'] + "alignment/minimap/from_filtering/{barcode}.minimap.sam",
        barcode=return_barcode_numbers(checkpoint_output))
def minimap_aligner_from_spoa(wildcards):
    #checkpoint_output = checkpoints.isONclustClusterFastq.get(**wildcards).output
    return config['results'] + "alignment/minimap/spoa.minimap.sam"
def vsearch_aligner(wildcards):
    isOnclustComplete = rules.isONclustClusterFastq.output.rule_complete
    return config['results'] + "alignment/vsearch/"
def id_reads(wildcards):
    checkpoint_output = rules.isONclustClusterFastq.output.rule_complete
    return [#config['results'] + "id_reads/id_reads.tsv",
            config['results'] + "id_reads/mapped_seq_id.csv",
            config['results'] + "id_reads/minimap_output.csv",
            config['results'] + "id_reads/mapped_consensus.csv"]
def spoa(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config['results'] + "spoa/consensus.sequences.fasta"
def nanoplot_basecall(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config['results'] + "visuals/nanoplot/basecall/"
def nanoplot_barcode_classified(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config['results'] + "visuals/nanoplot/barcode/classified"
def nanoplot_barcode_unclassified(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config['results'] + "visuals/nanoplot/barcode/unclassified"
def plotly_histogram_barcode(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config['results'] + "visuals/plotly/histograms/plotly.barcode.histogram.html"
def plotly_histogram_cutadapt(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config['results'] + "visuals/plotly/histograms/plotly.cutadapt.histogram.html"
def plotly_histogram_filtering(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config['results'] + "visuals/plotly/histograms/plotly.filtering.histogram.html"
def plotly_histogram_mapping(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config['results'] + "visuals/plotly/histograms/plotly.mapping.histogram.html"
def plotly_box_whisker(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config['results'] + "visuals/plotly/plotly.box.whisker.html"
FAST5_FILES = glob_wildcards(config['data'] + "fast5/" + "{fast5_file}.fast5").fast5_file

rule all:
    input:
        barcode_merge_files,                    # basecall, barcode, and merge files in a checkpoint
        nanoqc_basecall_data,                   # nanoqc after basecall
        nanoqc_barcode_classified,              # nanoqc classified barcodes
        nanoqc_barcode_unclassified,            # nanoqc unclassified barcodes
        cutadapt,                               # cutadapt on merged files
        filtering,                              # nanofilt on cutadapt files
        isONclust_pipeline,                     # cluster reads
        isONclust_cluster_fastq,                # clustering reads
        # guppy_aligner,                          # guppy
        minimap_aligner_from_filtering,         # minimap from filtering
        minimap_aligner_from_spoa,          # minimap from spoa clustering
        # vsearch_aligner,                       # vsearch
        id_reads,                              # id reads through python script
        IsoCon,                                # get consensus sequence
        spoa,                                  # partial order alignment
        nanoplot_basecall,                      # nanoplot for basecall files
        nanoplot_barcode_classified,            # nanoplot for classified barcode files
        nanoplot_barcode_unclassified,          # nanoplot for unclassified barcode files
        plotly_histogram_barcode,               # plotly barcode histogram
        plotly_histogram_cutadapt,              # plotly cutadapt histogram
        plotly_histogram_filtering,             # plotly filtering histogram
        plotly_histogram_mapping,               # plotly mapping histogram
        plotly_box_whisker,                     # plotly box and whisker plot

if config['basecall']:
    checkpoint basecall:
        input:
            config['data'] + "fast5"
        output:
            output = directory(config['results'] + "basecall/")
        params:
            config = os.environ['basecall_configuration'],
            callers = os.environ["num_basecallers"],
            threads_per_caller = os.environ["threads_per_caller"]
        shell:
            r"""
            echo Basecalling           
            
            guppy_basecaller \
            --config {params.config} \
            --input_path {input} \
            --save_path {output} \
            --num_callers {params.callers} \
            --cpu_threads_per_caller {params.threads_per_caller} \
            --recursive
            """

def barcode_input(wildcards):
    if config['basecall']:
        return rules.basecall.output[0]
    else:
        return config['data'] + "fastq"
checkpoint barcode:
    input:
        barcode_input
    output:
        output_directory = temp(directory(config['results'] + ".temp/barcodeTempOutput/")),
        barcode_complete_file = config['results'] + ".temp/barcodingDone"
    params:
        barcode_kit = os.environ['barcode_kit']
    shell:
        r"""
        guppy_barcoder \
        -i {input} \
        -s {output.output_directory} \
        --barcode_kits {params.barcode_kit} \
        --recursive 
        
        touch {output.barcode_complete_file}
        """

def merge_files_input(wildcards):
    return glob.glob(config['results'] + f".temp/barcodeTempOutput/{wildcards.barcode}/*.fastq")
rule merge_files:
    input:
        merge_files_input
    output:
        config['results'] + "barcode/{barcode}.merged.fastq"
    params:
        input_folder = config['results'] + ".barcodeTempOutput",
        save_folder = config['results'] + "barcode"
    shell:
        r"""
        for item in {input}; do
            cat $item >> {output}
        done
        """


def collate_basecall_fastq_files_input(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return glob.glob(config['results'] + "basecall/*.fastq")
def create_classified_unclassified_input(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return glob.glob(config['results'] + ".temp/barcodeTempOutput/**/*.fastq")
rule collate_basecall_fastq_files:
    input:
        collate_basecall_fastq_files_input
    output:
        temp(config['results'] + ".temp/basecall.temp.merged.files.fastq")
    shell:
        r"""        
        # concatenate each file in the params directory to the output file
        touch {output}
        for file in {input}; do
            cat "$file" >> {output}
        done
        """
rule create_classified_unclassified_file:
    input:
        create_classified_unclassified_input
    output:
        classified = temp(config['results'] + ".temp/barcode.classified.merged.temp.fastq"),
        unclassified = temp(config['results'] + ".temp/barcode.unclassified.merged.temp.fastq")
    shell:
        r"""
        for file in {input}; do
            if [[ $file == */barcode??/* ]]; then
                cat "$file" >> {output.classified}
            elif [[ $file == */unclassified/* ]]; then
                cat "$file" >> {output.unclassified}
            fi
        done
        """

rule NanoQCBasecall:
    input:
        rules.collate_basecall_fastq_files.output[0]
    output:
        directory(config['results'] + "visuals/nanoqc/basecall/")
    shell:
        r"""
        nanoQC -o {output} {input}
        """
rule NanoQCBarcode:
    input:
        classified = rules.create_classified_unclassified_file.output.classified,
        unclassified = rules.create_classified_unclassified_file.output.unclassified
    output:
        classified = directory(config['results'] + "visuals/nanoqc/barcode/classified"),
        unclassified = directory(config['results'] + "visuals/nanoqc/barcode/unclassified")
    shell:
        r"""
        nanoQC -o {output.classified} {input.classified}
        nanoQC -o {output.unclassified} {input.unclassified}
        """


rule NanoPlotBasecall:
    input:
        rules.collate_basecall_fastq_files.output[0]
    output:
        directory(config['results'] + "visuals/nanoplot/basecall/")
    run:
        # if we are not basecalling, the input file will have 0 lines
        # if this is the case, do nothing
        if len(open(str(input)).readlines()) != 0:
            subprocess.run(["NanoPlot", "--fastq", input, "-o", output])
        else:
            pass
rule NanoPlotBarcode:
    input:
        classified = rules.create_classified_unclassified_file.output.classified,
        unclassified = rules.create_classified_unclassified_file.output.unclassified

    output:
        classified = directory(config['results'] + "visuals/nanoplot/barcode/classified"),
        unclassified = directory(config['results'] + "visuals/nanoplot/barcode/unclassified")
    shell:
        r"""
        NanoPlot --fastq {input.classified} -o {output.classified}
        NanoPlot --fastq {input.unclassified} -o {output.unclassified}
        """


rule cutadapt:
    input:
        rules.merge_files.output[0]
    output:
        cutadapt_file = config['results'] + "cutadapt/{barcode}.cutadapt.fastq"
    params:
        three_prime_adapter = os.environ['cutadapt_trim_three_prime_adapter'],
        five_prime_adapter = os.environ['cutadapt_trim_five_prime_adapter'],
        error_rate = os.environ['cutadapt_trim_error_rate']

    shell:
        r"""
        cutadapt \
        --revcomp \
        --quiet \
        --cores 0 \
        --adapter {params.three_prime_adapter} \
        --front {params.five_prime_adapter} \
        --error-rate {params.error_rate} \
        --output {output} \
        {input}
        """
rule cutadaptDone:
    input:
        expand(rules.cutadapt.output, barcode=glob_wildcards(config['results'] + ".temp/barcodeTempOutput/{barcode}/*.fastq").barcode)
    output:
        touch(config['results'] + ".temp/cutadaptDone")
    shell:
        """
        # this is to ensure cutadapt is done before continuing.
        # Attempting to use cutadapt with a checkpoint results in the inability to fill the wildcard `barcode`
        """


rule filtering:
    input:
        rules.cutadapt.output[0]
    output:
        barcode_files = config['results'] + "filter/{barcode}.filter.fastq"
    params:
        min_length = os.environ['nanofilt_filtering_min'],
        max_length = os.environ['nanofilt_filtering_max'],
    shell:
        r"""
        touch {output}
        NanoFilt --length {params.min_length} --maxlength {params.max_length} {input} > {output}
        """


def merge_filtering_input(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    files = return_barcode_numbers(checkpoint_output)
    return expand(config['results'] + "filter/{barcode}.filter.fastq",
                  barcode=files)
rule merge_filtering_files:
    input:
        merge_filtering_input
    output:
        config['results'] + ".temp/merge.filtering.files.fastq"
    shell:
        r"""
        for file in {input}; do
            cat "$file" >> {output}
        done
        """



rule isOnClustPipeline:
    input:
        rules.merge_filtering_files.output[0]
    output:
        directory(config['results'] + "isONclust/pipeline/")
    shell:
        r"""
        # create a .tsv file
        
        isONclust --ont \
        --fastq {input} \
        --aligned_threshold 0.90 \
        --min_fraction 0.55 \
		--mapped_threshold 0.65 \
        --outfolder {output}
        """


rule isONclustClusterFastq:
    input:
        pipeline_output = rules.isOnClustPipeline.output[0],
        merged_filtering_reads = rules.merge_filtering_files.output[0]
    output:
        cluster_output = directory(config['results'] + "isONclust/cluster_fastq/"),
        rule_complete = config['results'] + ".temp/isONclustClusterFastqComplete"
    shell:
        r"""
        isONclust write_fastq --clusters {input.pipeline_output}/final_clusters.tsv \
        --fastq {input.merged_filtering_reads} \
        --outfolder {output.cluster_output} \
        --N 1
        
        touch {output.rule_complete}
        """



rule temp_spoa:
    input:
        complete_rule = rules.isONclustClusterFastq.output.rule_complete,
        cluster_data = expand(rules.isONclustClusterFastq.output.cluster_output + "{file}.fastq",
                              file=glob_wildcards(rules.isONclustClusterFastq.output.cluster_output + "{file}.fastq").file)
    output:
        temp_output = temp(directory(config['results'] + ".temp/spoa"))
    run:
        # make our temporary output folder
        os.mkdir(output.temp_output)

        for file in input.cluster_data:
            if ".fastq" in file:
                # convert fastq file names to fasta
                fastq_to_fasta = file[:-1] + "a"
                fastq_to_fasta = os.path.basename(fastq_to_fasta)

                # get our temporary output path generated and create the file
                temp_output = os.path.join(output.temp_output, fastq_to_fasta)

                # run spoa on the current isOnclust file, putting output into the temp spoa file
                with open(temp_output, 'w') as output_stream:
                    subprocess.run(["spoa", file, "-r", "0"], stdout=output_stream, universal_newlines=True)


rule spoa:
    input:
        rules.temp_spoa.output.temp_output
    output:
        config['results'] + "spoa/consensus.sequences.fasta"
    run:
        for file in os.listdir(str(input)):
            # we do not want the `.snakemake_timestamp` file to be included in this
            if ".snakemake_timestamp" not in file:
                file_path = os.path.join(str(input), file)

                # get file name without extension
                file_basename = file.split(".")[0]

                # overwrite first line in file with `>cluster_{file_basename}`
                # read lines
                file_lines = open(file_path, 'r').readlines()
                # replace first line
                file_lines[0] = f">cluster_{file_basename}\n"

                # append new lines to the output file
                open(str(output), 'a').writelines(file_lines)


rule guppy_aligner:
    input:
        rules.filtering.output[0]
    output:
        sam_files = touch(config['results'] + "alignment/guppy/sam_files/{barcode}.guppy.sam"),
        alignment_summary = touch(config['results'] + "alignment/guppy/alignment_summary/{barcode}.alignment.summary.csv"),
        log_file = touch(config['results'] + "alignment/guppy/logs/{barcode}.guppy.log"),
    params:
        barcode = "{barcode}",
        temp_dir = config['results'] + ".temp/guppy",
        alignment_reference = os.environ["alignment_path"]
    shell:
        r"""
        # move input files to our temp folder
        # this is required because guppy_aligner wants folders as input
        temp_input={params.temp_dir}/input/{params.barcode}
        temp_output={params.temp_dir}/output/{params.barcode}
        
        # in case any leftover barcode folders are present, we want to remove them, then recreate them
        rm -rf $temp_input; mkdir -p $temp_input
        rm -rf $temp_output        
        cp {input} $temp_input
        
        # run alignment
        guppy_aligner \
        --input_path $temp_input \
        --save_path $temp_output \
        --align_ref {params.alignment_reference} \
        --quiet 
        
        # move files from the temp output to the appropriate folders
        # a temp output allows for better organization
        mv $temp_output/{params.barcode}.filter.sam {output.sam_files}
        mv $temp_output/alignment_summary.txt {output.alignment_summary}
        mv $temp_output/read_processor_log*.log {output.log_file}
        
        # remove our temporary input and output files
        rm -rf $temp_input
        rm -rf $temp_output
        """



rule minimap_aligner_from_filtering:
    input:
        rules.filtering.output[0]
    output:
        config['results'] + "alignment/minimap/from_filtering/{barcode}.minimap.sam"
    params:
        alignment_reference = os.environ["alignment_path"]
    shell:
        r"""
        touch {output}
        minimap2 \
        -ax map-ont \
        {params.alignment_reference} \
        {input} > {output}
        """
rule minimap_aligner_from_spoa:
    input:
        rules.spoa.output[0]
    output:
        config['results'] + "alignment/minimap/spoa.minimap.sam"
    params:
        alignment_reference = os.environ["alignment_path"]
    shell:
        r"""
        touch {output}
        
        minimap2 \
        -ax map-ont \
        {params.alignment_reference} \
        {input} > {output}
        """



rule fq2fa:
    input:
        complete_rule = rules.isONclustClusterFastq.output.rule_complete,
        cluster_data = expand(rules.isONclustClusterFastq.output.cluster_output + "{file}.fastq",
                              file=glob_wildcards(rules.isONclustClusterFastq.output.cluster_output + "{file}.fastq").file)
    output:
        temp(directory(config['results'] + ".temp/vsearch/"))
    run:
        for item in input.cluster_data:
            file_name = os.path.basename(item)
            output_file = output[0] + file_name
            with open(output_file, 'w') as output_stream:
                subprocess.run(["seqkit", "fq2fa", item], stdout=output_stream, universal_newlines=True)
rule vsearch_aligner:
    input:
        rules.fq2fa.output
    output:
        directory(config['results'] + "alignment/vsearch/")
    params:
        alignment_reference = os.environ["alignment_path"]
    run:
        # {file}.fastq
        # vsearch.{file_number}.tsv
        for item in os.listdir(str(input)):
            # get just the file number
            file_number = os.path.basename(item)
            file_number = file_number.split(".")[0]

            # make our output path
            output_path = config['results'] + f"alignment/vsearch/vsearch.{file_number}.tsv"

            # call vsearch
            command = f"vsearch --sintax {input}{item} --tabbedout {output_path} --db {params.alignment_reference} --quiet"
            with open(output_path, 'w') as output_stream:
                subprocess.run(command.split(" "), stdout=output_stream, universal_newlines=True)


rule id_reads:
    input:
        filtering = expand(rules.filtering.output[0],
                           barcode=glob_wildcards(config['results'] + "filter/{barcode}.filter.fastq").barcode),
        clustering = rules.isOnClustPipeline.output[0],
        minimap = rules.minimap_aligner_from_spoa.output[0]
    output:
        #id_reads_tsc = config['results'] + "id_reads/id_reads.tsv",
        mapped_seq_id_csv = config['results'] + "id_reads/mapped_seq_id.csv",
        minimap_output = config['results'] + "id_reads/minimap_output.csv",
        mapped_consensus_csv = config['results'] + "id_reads/mapped_consensus.csv"

    params:
        results_folder = config['results']
    script:
        "scripts/id_reads.py"



rule IsoCon:
    input:
        merged_filter_files = rules.merge_filtering_files.output[0]
    output:
        directory(config['results'] + "isoCon/")
    shell:
        r"""
        IsoCon pipeline -fl_reads {input.merged_filter_files} -outfolder {output}
        """

def count_reads_barcode_input(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    files = return_barcode_numbers(checkpoint_output)
    merged_barcode_files = expand(rules.merge_files.output, barcode=files)
    return merged_barcode_files
def count_reads_cutadapt_input(wildcards):
    cutadapt_done = rules.cutadaptDone.output[0]
    cutadapt_output = rules.cutadapt.output[0]
    barcode_numbers = return_barcode_numbers(checkpoints.barcode.get(**wildcards).output[0])
    return expand(cutadapt_output, barcode=barcode_numbers)
def count_reads_filtering_input(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    files = return_barcode_numbers(checkpoint_output)
    return expand(config['results'] + "filter/{barcode}.filter.fastq", barcode=files)
def count_minimap_reads(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    files = return_barcode_numbers(checkpoint_output)
    return expand(config['results'] + "alignment/minimap/from_filtering/{barcode}.minimap.sam", barcode=files)
rule count_reads_barcode:
    input:
        count_reads_barcode_input
    output:
        config['results'] + "count_reads/count.reads.barcode.csv"
    params:
        process = "barcode"
    script:
        "scripts/CountReads.py"
rule count_reads_cutadapt:
    input:
        count_reads_cutadapt_input
    output:
        config['results'] + "count_reads/count.reads.cutadapt.csv"
    params:
        process = "cutadapt"
    script:
        "scripts/CountReads.py"
rule count_reads_filtering:
    input:
        count_reads_filtering_input
    output:
        config['results'] + "count_reads/count.reads.filter.csv"
    params:
        process = "filtering"
    script:
        "scripts/CountReads.py"
rule count_reads_mapping:
    input:
        count_minimap_reads
    output:
        config['results'] + "count_reads/count.reads.mapping.csv"
    params:
        process = "mapping"
    script:
        "scripts/CountReads.py"



rule plotly_barcode_histogram:
    input:
        rules.count_reads_barcode.output[0]
    output:
        config['results'] + "visuals/plotly/histograms/plotly.barcode.histogram.html"
    params:
        sub_title = "Performed after Merging Files"
    script:
        "scripts/PlotlyHistogram.py"
rule plotly_cutadapt_histogram:
    input:
        rules.count_reads_cutadapt.output[0]
    output:
        config['results'] + "visuals/plotly/histograms/plotly.cutadapt.histogram.html"
    params:
        sub_title = "Performed after Cutadapt",
    script:
        "scripts/PlotlyHistogram.py"
rule plotly_filtering_histogram:
    input:
        rules.count_reads_filtering.output[0]
    output:
        config['results'] + "visuals/plotly/histograms/plotly.filtering.histogram.html"
    params:
        sub_title = "Performed after Filtering",
    script:
        "scripts/PlotlyHistogram.py"
rule plotly_mapping_histogram:
    input:
        rules.count_reads_mapping.output[0]
    output:
        config['results'] + "visuals/plotly/histograms/plotly.mapping.histogram.html"
    params:
        sub_title = "Performed after Mapping"
    script:
        "scripts/PlotlyHistogram.py"
rule plotly_box_whisker_generation:
    input:
        rules.count_reads_barcode.output[0],
        rules.count_reads_cutadapt.output[0],
        rules.count_reads_filtering.output[0],
        rules.count_reads_mapping.output[0]
    output:
        config['results'] + "visuals/plotly/plotly.box.whisker.html"
    script:
        "scripts/PlotlyBoxWhisker.py"

