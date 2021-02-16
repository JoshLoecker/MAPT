import os
from pathlib import Path
import glob
from multiprocessing import cpu_count
from pprint import pprint
import shutil
import pandas as pd
import numpy as np
import gzip

configfile: "config.yaml"

def return_barcode_numbers(path: str):
    """
    This function will return a list of barcode numbers from the input path
    The input path should contain barcode folder names
    Returns: A list of strings containing directory names
    """

    barcode_set = set()
    for folder in os.scandir(path):
        if folder.is_dir():
            barcode_set.add(folder.name)
    return barcode_set
def barcode_merge_files(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    barcodes = set()  # a set is like a list but only stores unique values
    for folder in os.scandir(checkpoint_output):
        if folder.is_dir():
            barcodes.add(folder.name)

    return_merged_files = [config["results"] + "barcode/" + barcode + ".merged.fastq.gz" for barcode in barcodes]
    return return_merged_files
def nanoqc_basecall_data(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config["results"] + "visuals/nanoqc/basecall/"
def nanoqc_barcode_classified(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config["results"] + "visuals/nanoqc/barcode/classified"
def nanoqc_barcode_unclassified(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config["results"] + "visuals/nanoqc/barcode/unclassified"
def cutadapt(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return expand(
        config["results"] + "cutadapt/{barcode}.cutadapt.fastq.gz",
        barcode=return_barcode_numbers(config["results"] + ".temp/barcodeTempOutput"))
def filtering(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return expand(config["results"] + "filter/{barcode}.filter.fastq.gz",
        barcode=return_barcode_numbers(config["results"] + ".temp/barcodeTempOutput"))
def isONclust_pipeline(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config["results"] + "isONclust/pipeline/"
def isONclust_cluster_fastq(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config["results"] + "isONclust/cluster_fastq/"
def minimap_aligner_from_filtering(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return expand(
        config["results"] + "alignment/minimap/from_filtering/{barcode}.minimap.sam",
        barcode=return_barcode_numbers(config["results"] + ".temp/barcodeTempOutput"))
def minimap_aligner_from_spoa(wildcards):
    return config["results"] + "alignment/minimap/from_spoa/spoa.minimap.sam"
def id_reads(wildcards):
    checkpoint_output = rules.isONclustClusterFastq.output.rule_complete
    return [config["results"] + "id_reads/mapped_reads/mapped_seq_id.csv",
            config["results"] + "id_reads/mapped_reads/minimap_output.csv",
            config["results"] + "id_reads/mapped_reads/mapped_consensus.csv"]
def spoa(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config["results"] + "spoa/consensus.sequences.fasta.gz"
def nanoplot_basecall(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config["results"] + "visuals/nanoplot/basecall/"
def nanoplot_barcode_classified(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config["results"] + "visuals/nanoplot/barcode/classified"
def nanoplot_barcode_unclassified(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config["results"] + "visuals/nanoplot/barcode/unclassified"
def plotly_histogram_barcode(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config["results"] + "visuals/plotly/histograms/plotly.barcode.histogram.html"
def plotly_histogram_cutadapt(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config["results"] + "visuals/plotly/histograms/plotly.cutadapt.histogram.html"
def plotly_histogram_filtering(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config["results"] + "visuals/plotly/histograms/plotly.filtering.histogram.html"
def plotly_histogram_mapping(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config["results"] + "visuals/plotly/histograms/plotly.mapping.histogram.html"
def plotly_box_whisker(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    return config["results"] + "visuals/plotly/plotly.box.whisker.html"

FAST5_FILES = glob_wildcards(config["basecall_files"] + "{fast5_file}.fast5").fast5_file


rule all:
    input:
        config["results"] + ".temp/completeRules/basecallComplete",# basecalling
        barcode_merge_files, # barcode and merge files
        cutadapt, # cutadapt on merged files
        filtering, # nanofilt on cutadapt files
        nanoqc_basecall_data,# nanoqc after basecall
        nanoqc_barcode_classified,# nanoqc classified barcodes
        nanoqc_barcode_unclassified,# nanoqc unclassified barcodes
        isONclust_pipeline,# cluster reads
        isONclust_cluster_fastq,# clustering reads
        minimap_aligner_from_filtering, # minimap from filtering
        minimap_aligner_from_spoa, # minimap from spoa clustering
        id_reads, # id reads through python script
        config["results"] + "id_reads/filter_id_reads/withinDivergence.csv",  # filter id_reads that are inside divergence
        config["results"] + "id_reads/filter_id_reads/outsideDivergence.csv",  # filter id_reads that are outside divergence
        config["results"] + "id_reads/filter_id_reads/nanDivergence.csv",  # filter id_reads that have no divergence
        config["results"] + "id_reads/OTU/withinDivergenceOTU.csv",  # create OTU table of values within divergence
        config["results"] + "id_reads/OTU/outsideDivergenceOTU.csv",  # create OTU table of values outside divergence
        config["results"] + "id_reads/OTU/nanDivergenceOTU.csv",  # create OTU table of id_reads that have no divergence
        config["results"] + ".temp/completeRules/RemoveLowClustersComplete",  # remove clusters with low reads
        config["results"] + "id_reads/simple_mapped_reads/simpleMappedWithinDivergence.csv",  # create simplified mapped_seq_id csv within divergence value
        config["results"] + "id_reads/simple_mapped_reads/simpleMappedOutsideDivergence.csv",  # create simplified mapped_seq_id csv outside divergence value
        config["results"] + "id_reads/simple_mapped_reads/simpleMappedNaNDivergence.csv",  # create simplified mapped_seq_id csv without divergence value
        config["results"] + "id_reads/cluster_summary/clusterSummaryWithinDivergence.csv",  # create a cluster summary file; csv files with data within/outside/no divergence data are created
        config["results"] + "id_reads/cluster_summary/clusterSummaryOutsideDivergence.csv",
        config["results"] + "id_reads/cluster_summary/clusterSummaryNaNDivergence.csv",
        spoa, # partial order alignment
        nanoplot_basecall,# nanoplot for basecall files
        nanoplot_barcode_classified,# nanoplot for classified barcode files
        nanoplot_barcode_unclassified,# nanoplot for unclassified barcode files
        plotly_histogram_barcode,# plotly barcode histogram
        plotly_histogram_cutadapt,# plotly cutadapt histogram
        plotly_histogram_filtering,# plotly filtering histogram
        plotly_histogram_mapping,# plotly mapping histogram
        plotly_box_whisker  # plotly box and whisker plot



"""
We are placing the output in the parameters section instead of the output section because snakemake will delete the output folder
    if the job is terminated (i.e. out of time on SciNet)
This is the only available option of not deleting output when snakemake is terminated
From: https://stackoverflow.com/questions/55419603/how-to-prevent-snakemake-from-deleting-output-folders-from-failed-jobs
"""
if config["basecall"]["perform_basecall"]:
    checkpoint basecall:
        input:
            config["basecall_files"]
        output:
            rule_complete=config["results"] + ".temp/completeRules/basecallComplete"
        params:
            guppy_container=config["guppy_container"],
            config=config["basecall"]["configuration"],
            output=config["results"] + "basecall"
        shell:
            r"""
            # Try to resume guppy_basecaller, otherwise simply execute guppy_basecaller
            
            command="singularity exec --nv {params.guppy_container}  \
                guppy_basecaller \
                --config {params.config} \
                --input_path {input} \
                --save_path {params.output} \
                --device 'cuda:all' \
                --recursive"

            # try to resume basecalling
            # if this fails, perform normal basecalling
            # this is ugly, but I am not sure of any other method to try resuming basecalling, and trying again without resuming
            eval "$command --resume || $command"
            
            # zip output
            gzip --best -f {params.output}/*.fastq
            touch {output.rule_complete}
            """

    def collate_basecall_fastq_files_input(wildcards):
        checkpoint_output = checkpoints.basecall.get(**wildcards).output[0]
        basecall_output = config["results"] + "basecall"
        return glob.glob(basecall_output + "/*.fastq.gz")
    rule collate_basecall_fastq_files:
        input:
            collate_basecall_fastq_files_input
        output:
            fastq_gz = temp(config["results"] + ".temp/basecall.merged.files.fastq.gz")
        params:
            temp_fastq = config["results"] + ".temp/basecall.merged.files.fastq"
        shell:
            r"""
            echo HERE
            echo {input}
            for file in {input}; do
                gunzip --stdout "$file" >> {params.temp_fastq}
            done
            
            gzip --best -f {params.temp_fastq}
            """


def barcode_input(wildcards):
    if config["basecall"]["perform_basecall"]:
        output = checkpoints.basecall.get(**wildcards).output
        return output
    else:
        return config["barcode_files"]
checkpoint barcode:
    input:
        barcode_input
    output:
        output=directory(config["results"] + ".temp/barcodeTempOutput/"),
        barcode_complete_file=config["results"] + ".temp/completeRules/barcodingComplete"
    params:
        guppy_container=config["guppy_container"],
        barcode_kit=config["barcode"]["kit"],
        basecall_output=config["results"] + "basecall/"
    shell:
        r"""
        command="singularity exec --nv {params.guppy_container} \
        guppy_barcoder \
        --input_path {params.basecall_output} \
        --save_path {output.output} \
        --barcode_kits {params.barcode_kit} \
        --recursive"
        
        local_command="guppy_barcoder \
        --input_path {params.basecall_output} \
        --save_path {output.output} \
        --barcode_kits {params.barcode_kit} \
        --recursive"
        
        # try to execute barcoding with a GPU, otherwise continue without
        eval "$command --device 'cuda:all' || $command || $local_command"
        
        touch {output.barcode_complete_file}
        """


def merge_files_input(wildcards):
    return glob.glob(config["results"] + f".temp/barcodeTempOutput/{wildcards.barcode}/*.fastq")
rule merge_files:
    input:
        merge_files_input
    output:
        config["results"] + "barcode/{barcode}.merged.fastq.gz"
    params:
        temp_file=config["results"] + "barcode/{barcode}.merged.fastq",
        barcode_output=config["results"] + ".temp/barcodeTempOutput"
    shell:
        r"""
        for item in {input}; do
            cat $item >> {params.temp_file}
            gzip --best -f {params.temp_file}
        done
        """
rule merge_files_done:
    input:
        expand(rules.merge_files.output,
            barcode=glob_wildcards(config["results"] + ".temp/barcodeTempOutput/{barcode}/*.fastq").barcode)
    output:
        touch(config["results"] + ".temp/completeRules/mergeFilesDone")
    shell:
        """
        # this is to ensure filtering is done before continuing.
        # Attempting to use filtering with a checkpoint results in the inability to fill the wildcard `barcode`
        """


def merge_barcode_files(wildcards):
    checkpoint_output = checkpoints.barcode.get(**wildcards).output
    barcode_numbers = return_barcode_numbers(config["results"] + ".temp/barcodeTempOutput/")
    return expand(config["results"] + "barcode/{barcode}.merged.fastq.gz", barcode=barcode_numbers)
rule create_classified_unclassified_file:
    input:
        merge_barcode_files
    output:
        classified=temp(config["results"] + ".temp/barcode.classified.merged.fastq.gz"),
        unclassified=temp(config["results"] + ".temp/barcode.unclassified.merged.fastq.gz")
    params:
        classified = temp(config["results"] + ".temp/barcode.classified.merged.fastq"),
        unclassified = temp(config["results"] + ".temp/barcode.unclassified.merged.fastq")
    shell:
        r"""
        # The '=~' operator is a regex match
        for file in {input}; do
            if [[ "$file" =~ barcode[0-9]{{2}}.merged ]]; then
                gunzip --stdout >> {params.classified}
            elif [[ "$file" =~ unclassified.merged ]]; then
                gunzip --stdout >> {params.unclassified}
            fi
        done
        
        gzip --best -f {params.classified}
        gzip --best -f {params.unclassified}
        
        """

rule NanoQCBasecall:
    input:
        rules.collate_basecall_fastq_files.output.fastq_gz
    output:
        directory(config["results"] + "visuals/nanoqc/basecall/")
    shell:
        r"""
        nanoQC -o {output} {input}
        """

rule NanoQCBarcode:
    input:
        classified=rules.create_classified_unclassified_file.output.classified,
        unclassified=rules.create_classified_unclassified_file.output.unclassified
    output:
        classified=directory(config["results"] + "visuals/nanoqc/barcode/classified"),
        unclassified=directory(config["results"] + "visuals/nanoqc/barcode/unclassified")
    shell:
        r"""
        nanoQC -o {output.classified} {input.classified}
        nanoQC -o {output.unclassified} {input.unclassified}
        """

if config["basecall"]["perform_basecall"]:
    rule NanoPlotBasecall:
        input:
            rules.collate_basecall_fastq_files.output.fastq_gz
        output:
            directory(config["results"] + "visuals/nanoplot/basecall/")
        params:
            temp_gunzip_file=config["results"] + ".temp/nanoplot.basecall.gunzip"
        shell:
            r"""
            NanoPlot --fastq {input} -o {output}
            # gunzip --stdout {input} > {params.temp_gunzip_file}
            # NanoPlot --fastq {params.temp_gunzip_file} -o {output}
            # rm {params.temp_gunzip_file}
            """

rule NanoPlotBarcode:
    input:
        classified=rules.create_classified_unclassified_file.output.classified,
        unclassified=rules.create_classified_unclassified_file.output.unclassified

    output:
        classified=directory(config["results"] + "visuals/nanoplot/barcode/classified"),
        unclassified=directory(config["results"] + "visuals/nanoplot/barcode/unclassified")
    shell:
        r"""
        NanoPlot --fastq {input.classified} -o {output.classified}
        NanoPlot --fastq {input.unclassified} -o {output.unclassified}
        """


rule cutadapt:
    input:
        rules.merge_files.output[0]
    output:
        cutadapt_file=config["results"] + "cutadapt/{barcode}.cutadapt.fastq.gz"
    params:
        three_prime_adapter=config["cutadapt"]["three_prime_adapter"],
        five_prime_adapter=config["cutadapt"]["five_prime_adapter"],
        error_rate=config["cutadapt"]["error_rate"],
        temp_fastq=config["results"] + "cutadapt/{barcode}.cutadapt.fastq"

    shell:
        r"""
        cutadapt \
        --revcomp \
        --quiet \
        --cores 0 \
        --adapter {params.three_prime_adapter} \
        --front {params.five_prime_adapter} \
        --error-rate {params.error_rate} \
        --output {params.temp_fastq} \
        {input}
        
        gzip --best -f {params.temp_fastq}
        """
rule cutadaptDone:
    input:
        expand(rules.cutadapt.output,
            barcode=glob_wildcards(config["results"] + ".temp/barcodeTempOutput/{barcode}/*.fastq").barcode)
    output:
        touch(config["results"] + ".temp/completeRules/cutadaptDone")
    shell:
        """
        # this is to ensure cutadapt is done before continuing.
        # Attempting to use cutadapt with a checkpoint results in the inability to fill the wildcard `barcode`
        """


rule filtering:
    input:
        rules.cutadapt.output.cutadapt_file
    output:
        filtering_files=config["results"] + "filter/{barcode}.filter.fastq.gz"
    params:
        min_quality=config["nanofilt"]["min_quality"],
        min_length=config["nanofilt"]["min_filter"],
        max_length=config["nanofilt"]["max_filter"],
        temp_fastq=config["results"] + "filter/{barcode}.filter.fastq"
    shell:
        r"""
        touch {output}
        
        gunzip --stdout {input} > {params.temp_fastq}
        
        NanoFilt \
        --quality {params.min_quality} \
        --length {params.min_length} \
        --maxlength {params.max_length} {params.temp_fastq} > {output.filtering_files}
        
        gzip --best -f {params.temp_fastq}
        """

def merge_filtering_input(wildcards):
    checkpoint_output = checkpoints.barcode.get(**wildcards).output
    barcode_numbers = return_barcode_numbers(config["results"] + ".temp/barcodeTempOutput/")
    # return expand(config["results"] + "filter/{barcode}.filter.fastq", barcode=barcode_numbers)
    return expand(rules.filtering.output.filtering_files, barcode=barcode_numbers)
rule merge_filtering:
    input:
        merge_filtering_input
    output:
        merged_filtering = config["results"] + ".temp/merge.filtering.files.fastq.gz"
    params:
        merged_filtering = config["results"] + ".temp/merge.filtering.files.fastq"
    shell:
        r"""
        for file in {input}; do
            cat "$file" >> {params.merged_filtering}
        done
        gzip --best -f {params.merged_filtering}
        """

rule isONClustPipeline:
    input:
        rules.merge_filtering.output.merged_filtering
    output:
        data=directory(config["results"] + "isONclust/pipeline/"),
        rule_complete=config["results"] + ".temp/completeRules/isONClustPipelineComplete"
    params:
        aligned_threshold=config["isONclust"]["aligned_threshold"],
        min_fraction=config["isONclust"]["min_fraction"],
        mapped_threshold=config["isONclust"]["mapped_threshold"],
        temp_fastq_input = config["results"] + ".temp/merge.filtering.files.isONclustPipeline.fastq"
    shell:
        r"""
        gunzip --stdout {input} > {params.temp_fastq_input}
        
        isONclust --ont \
        --fastq {params.temp_fastq_input} \
        --aligned_threshold {params.aligned_threshold} \
        --min_fraction {params.min_fraction} \
        --mapped_threshold {params.mapped_threshold} \
        --outfolder {output.data}
        
        rm -rf {params.temp_fastq_input}
        touch {output.rule_complete}
        """

rule isONclustClusterFastq:
    input:
        pipeline_output=rules.isONClustPipeline.output.data,
        isONClustComplete=rules.isONClustPipeline.output.rule_complete,
        merged_filtering_reads=rules.merge_filtering.output.merged_filtering
    output:
        cluster_output=directory(config["results"] + "isONclust/cluster_fastq/"),
        rule_complete=config["results"] + ".temp/completeRules/isONclustClusterFastqComplete"
    params:
        min_quality=config["nanofilt"]["min_quality"],  # use same quality as NanoFilt (i.e. rule filtering)
        temp_fastq_input = config["results"] + ".temp/merge.filtering.files.isONclustClusterFastq.fastq"
    shell:
        r"""
        gunzip --stdout {input.merged_filtering_reads} > {params.temp_fastq_input}
        
        isONclust \
        --q "{params.min_quality}" \
        write_fastq \
        --fastq {params.temp_fastq_input} \
        --outfolder "{output.cluster_output}" \
        --clusters "{input.pipeline_output}/final_clusters.tsv"
                
        gzip --best -f {output.cluster_output}/*.fastq
        rm -f {params.temp_fastq_input}
        touch {output.rule_complete}
        """


rule remove_low_reads:
    input:
        previous_rule_complete=rules.isONclustClusterFastq.output.rule_complete,
        cluster_data=expand(rules.isONclustClusterFastq.output.cluster_output + "{file}.fastq.gz",
            file=glob_wildcards(rules.isONclustClusterFastq.output.cluster_output + "{file}.fastq.gz").file)
    output:
        rule_complete=touch(config["results"] + ".temp/completeRules/RemoveLowClustersComplete")
    params:
        cluster_cutoff=config["cluster"]["min_reads_per_cluster"],
        output_folder=config["results"] + "TooFewReadsInCluster",
        temp_fastq = config["results"] + ".temp/remove_low_reads"
    run:
        os.makedirs(params.output_folder,exist_ok=True)
        os.makedirs(params.temp_fastq, exist_ok=True)
        input_data = str(input.cluster_data).split(" ")
        for file in input_data:
            # unzip input file
            temp_file_name = "{params.temp_fastq}/{file}.fastq"
            shell(f"gunzip --stdout {input.cluster_data} > {temp_file_name}")
            lines_in_file = open(temp_file_name, "r").readlines()
            count_lines_in_file = len(lines_in_file)
            count_reads_in_file = count_lines_in_file / 4

            # keep files that are equal to OR GREATER than the cutoff
            if count_reads_in_file < params.cluster_cutoff:
                # copy the file, then remove it
                # doing this prevents errors about the destination file already existing
                shutil.copy(src=file,dst=str(params.output_folder))
                os.remove(file)

            # remove unzipped file
            shell("rm {temp_file_name}")


rule temp_spoa:
    input:
        complete_rule=rules.remove_low_reads.output.rule_complete,
        cluster_data=expand(rules.isONclustClusterFastq.output.cluster_output + "{file}.fastq.gz",
            file=glob_wildcards(rules.isONclustClusterFastq.output.cluster_output + "{file}.fastq.gz").file)
    output:
        temp_output=directory(config["results"] + ".temp/spoa"),
        rule_complete=config["results"] + ".temp/completeRules/tempSpoaComplete"
    shell:
        r"""
        mkdir -p {output.temp_output}
        for file in {input.cluster_data}; do
            file_basename=$(basename -- "$file")
            output_file="{output.temp_output}/$file_basename"            
            spoa "$file" -r 0 > "$output_file"
        done
        
        touch {output.rule_complete}
        """

rule spoa:
    input:
        rule_complete=rules.temp_spoa.output.rule_complete,
        data=rules.temp_spoa.output.temp_output
    output:
        config["results"] + "spoa/consensus.sequences.fasta.gz"
    params:
        temp_fasta=config["results"] + "spoa/consensus.sequences.fasta"
    run:
        for file in os.listdir(str(input.data)):
            # we do not want the `.snakemake_timestamp` file to be included in this
            if ".snakemake_timestamp" not in file:
                file_path = os.path.join(str(input.data),file)

                # get file name without extension
                cluster_number = file.split(".")[0]

                # overwrite first line in file with `>cluster_{cluster_number}`
                file_lines = open(file_path,"r").readlines()
                file_lines[0] = f">cluster_{cluster_number}\n"

                # append new lines to output file
                open(str(params.temp_fasta),"a").writelines(file_lines)

                # remove file under results/.temp/spoa/*.fastq, it is no longer required
                os.remove(file_path)

        shell("gzip --best -f {params.temp_fasta}")


rule minimap_aligner_from_filtering:
    input:
        rules.filtering.output[0]
    output:
        config["results"] + "alignment/minimap/from_filtering/{barcode}.minimap.sam"
    params:
        alignment_reference=config["reference_database"],
        temp_fastq = config["results"] + ".temp/minimap/{barcode}.temp.minimap.fastq"
    shell:
        r"""
        mkdir -p {config[results]}.temp/minimap
        touch {output}
        touch {params.temp_fastq}
        
        gunzip --stdout {input} > {params.temp_fastq}
        
        minimap2 \
        -ax map-ont \
        {params.alignment_reference} \
        {params.temp_fastq} > {output}
        
        rm -f {params.temp_fastq}
        """

rule minimap_aligner_from_spoa:
    input:
        rules.spoa.output[0]
    output:
        config["results"] + "alignment/minimap/from_spoa/spoa.minimap.sam"
    params:
        alignment_reference=config["reference_database"],
        temp_fastq = config["results"] + ".temp/minimap/temp_spoa.fastq"
    shell:
        r"""
        touch {output}
        
        gunzip --stdout {input} > {params.temp_fastq}

        minimap2 \
        -ax map-ont \
        {params.alignment_reference} \
        {params.temp_fastq} > {output}
        
        rm -r {params.temp_fastq}
        """


rule id_reads:
    input:
        filtering=expand(rules.filtering.output[0],
            barcode=glob_wildcards(config["results"] + "filter/{barcode}.filter.fastq").barcode),
        clustering=rules.isONClustPipeline.output[0],
        minimap=rules.minimap_aligner_from_spoa.output[0]
    output:
        #id_reads_tsc = config["results"] + "id_reads/id_reads.tsv",
        mapped_seq_id_csv=config["results"] + "id_reads/mapped_reads/mapped_seq_id.csv",
        minimap_output=config["results"] + "id_reads/mapped_reads/minimap_output.csv",
        mapped_consensus_csv=config["results"] + "id_reads/mapped_reads/mapped_consensus.csv"

    params:
        results_folder=config["results"]
    script:
        "scripts/id_reads.py"


rule filter_id_reads_mapped_sequence:
    input:
        csv=rules.id_reads.output.mapped_seq_id_csv
    output:
        within_divergence=config["results"] + "id_reads/filter_id_reads/withinDivergence.csv",
        outside_divergence=config["results"] + "id_reads/filter_id_reads/outsideDivergence.csv",
        nan_divergence=config["results"] + "id_reads/filter_id_reads/nanDivergence.csv"
    params:
        divergence_threshold=config["cluster"]["divergence_threshold"]
    run:
        data_frame = pd.read_csv(input.csv,delimiter=",",header=0)
        header_data = data_frame.columns.values

        # create three pandas dataframes. One within bounds, one outside bounds, and one with NaN data
        within_bounds_data = data_frame.loc[data_frame["divergence"] <= params.divergence_threshold]
        outside_bounds_data = data_frame.loc[data_frame["divergence"] > params.divergence_threshold]
        nan_data = data_frame.loc[pd.isnull(data_frame["divergence"])]

        # write data to respective csv
        within_bounds_data.to_csv(path_or_buf=str(output.within_divergence),header=header_data,index=False)
        outside_bounds_data.to_csv(path_or_buf=str(output.outside_divergence),header=header_data,index=False)
        nan_data.to_csv(path_or_buf=str(output.nan_divergence),header=header_data,index=False)


rule otu_from_filter_id_reads:
    input:
        within_divergence=rules.filter_id_reads_mapped_sequence.output.within_divergence,
        outside_divergence=rules.filter_id_reads_mapped_sequence.output.outside_divergence,
        nan_divergence=rules.filter_id_reads_mapped_sequence.output.nan_divergence
    output:
        within_divergence_otu=config["results"] + "id_reads/OTU/withinDivergenceOTU.csv",
        outside_divergence_otu=config["results"] + "id_reads/OTU/outsideDivergenceOTU.csv",
        nan_divergence_otu=config["results"] + "id_reads/OTU/nanDivergenceOTU.csv"
    script:
        "scripts/generateOTU.py"


rule make_simple_mapped_sequence_id:
    input:
        rules.id_reads.output.mapped_seq_id_csv
    output:
        within_divergence=config["results"] + "id_reads/simple_mapped_reads/simpleMappedWithinDivergence.csv",
        outside_divergence=config["results"] + "id_reads/simple_mapped_reads/simpleMappedOutsideDivergence.csv",
        nan_divergence=config["results"] + "id_reads/simple_mapped_reads/simpleMappedNaNDivergence.csv"
    params:
        divergence_threshold=config["cluster"]["divergence_threshold"]
    script:
        "scripts/simpleMappedSequenceID.py"


rule cluster_summary:
    input:
        rules.id_reads.output.mapped_seq_id_csv
    output:
        within_divergence=config["results"] + "id_reads/cluster_summary/clusterSummaryWithinDivergence.csv",
        outside_divergence=config["results"] + "id_reads/cluster_summary/clusterSummaryOutsideDivergence.csv",
        nan_divergence=config["results"] + "id_reads/cluster_summary/clusterSummaryNaNDivergence.csv"
    params:
        divergence_threshold=config["cluster"]["divergence_threshold"]
    script:
        "scripts/clusterSummary.py"


def count_reads_barcode_input(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    files = return_barcode_numbers(config["results"] + ".temp/barcodeTempOutput")
    merged_barcode_files = expand(rules.merge_files.output,barcode=files)
    return merged_barcode_files


def count_reads_cutadapt_input(wildcards):
    cutadapt_done = rules.cutadaptDone.output[0]
    cutadapt_output = rules.cutadapt.output[0]
    barcode_numbers = return_barcode_numbers(config["results"] + ".temp/barcodeTempOutput")
    return expand(cutadapt_output,barcode=barcode_numbers)


def count_reads_filtering_input(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    files = return_barcode_numbers(config["results"] + ".temp/barcodeTempOutput")
    return expand(config["results"] + "filter/{barcode}.filter.fastq",barcode=files)


def count_minimap_reads(wildcards):
    barcode_done = checkpoints.barcode.get(**wildcards).output[1]
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    files = return_barcode_numbers(config["results"] + ".temp/barcodeTempOutput")
    return expand(config["results"] + "alignment/minimap/from_filtering/{barcode}.minimap.sam",barcode=files)


rule count_reads_barcode:
    input:
        count_reads_barcode_input
    output:
        config["results"] + "count_reads/count.reads.barcode.csv"
    params:
        process="barcode"
    script:
        "scripts/CountReads.py"

rule count_reads_cutadapt:
    input:
        count_reads_cutadapt_input
    output:
        config["results"] + "count_reads/count.reads.cutadapt.csv"
    params:
        process="cutadapt"
    script:
        "scripts/CountReads.py"

rule count_reads_filtering:
    input:
        count_reads_filtering_input
    output:
        config["results"] + "count_reads/count.reads.filter.csv"
    params:
        process="filtering"
    script:
        "scripts/CountReads.py"

rule count_reads_mapping:
    input:
        count_minimap_reads
    output:
        config["results"] + "count_reads/count.reads.mapping.csv"
    params:
        process="mapping"
    script:
        "scripts/CountReads.py"


rule plotly_barcode_histogram:
    input:
        rules.count_reads_barcode.output[0]
    output:
        config["results"] + "visuals/plotly/histograms/plotly.barcode.histogram.html"
    params:
        sub_title="Performed after Merging Files"
    script:
        "scripts/PlotlyHistogram.py"

rule plotly_cutadapt_histogram:
    input:
        rules.count_reads_cutadapt.output[0]
    output:
        config["results"] + "visuals/plotly/histograms/plotly.cutadapt.histogram.html"
    params:
        sub_title="Performed after Cutadapt"
    script:
        "scripts/PlotlyHistogram.py"

rule plotly_filtering_histogram:
    input:
        rules.count_reads_filtering.output[0]
    output:
        config["results"] + "visuals/plotly/histograms/plotly.filtering.histogram.html"
    params:
        sub_title="Performed after Filtering"
    script:
        "scripts/PlotlyHistogram.py"

rule plotly_mapping_histogram:
    input:
        rules.count_reads_mapping.output[0]
    output:
        config["results"] + "visuals/plotly/histograms/plotly.mapping.histogram.html"
    params:
        sub_title="Performed after Mapping"
    script:
        "scripts/PlotlyHistogram.py"

rule plotly_box_whisker_generation:
    input:
        rules.count_reads_barcode.output[0],
        rules.count_reads_cutadapt.output[0],
        rules.count_reads_filtering.output[0],
        rules.count_reads_mapping.output[0]
    output:
        config["results"] + "visuals/plotly/plotly.box.whisker.html"
    script:
        "scripts/PlotlyBoxWhisker.py"
