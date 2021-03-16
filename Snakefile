import glob
import random
import os
import re
import shutil
from pprint import pprint
import pandas as pd
# TODO: gzip all fastq files
configfile: "config.yaml"

def return_barcodes(wildcards):
    """
    This will return barcode numbers from the barcode output folder
    Simple pass through the `wildcards` parameter
    :param wildcards: The wildcards parameter used in checkpoint outputs
    :return:
    """
    checkpoint_output = checkpoints.barcode.get(**wildcards).output
    barcodes = set()
    for folder in os.scandir(os.path.join(config["results"], ".temp/barcode")):
        if folder.is_dir():
            barcodes.add(folder.name)
    return barcodes
def merge_barcodes(wildcards):
    barcodes = return_barcodes(wildcards)
    return expand(os.path.join(config["results"], "barcode/{barcode}.merged.fastq"), barcode=barcodes)
def basecall_visuals(wildcards):
    """
    If we are running basecalling, then NanoQC and NanoPlot rules should be ran
    :return:
    """
    if config["basecall"]["perform_basecall"]:
        return [os.path.join(config["results"], "visuals/nanoqc/basecall/"),
                os.path.join(config["results"], "visuals/nanoplot/basecall/")]
    else:
        return []
def minimap_from_filter(wildcards):
    barcodes = return_barcodes(wildcards)
    return expand(os.path.join(config["results"], "alignment/minimap/from_filtering/{barcode}.minimap.sam"), barcode=barcodes)


rule all:
    input:
        os.path.join(config["results"], "visuals/nanoplot/barcode/classified"),
        os.path.join(config["results"], "visuals/nanoplot/barcode/unclassified"),

        os.path.join(config["results"], "visuals/nanoqc/barcode/classified"),
        os.path.join(config["results"], "visuals/nanoqc/barcode/unclassified"),

        basecall_visuals,
        os.path.join(config["results"], "isONclust/cluster_fastq"),
        os.path.join(config["results"],"isONclust/pipeline"),
        os.path.join(config["results"], "LowClusterReads"),
        os.path.join(config["results"], "spoa/consensus.sequences.fasta"),
        minimap_from_filter,
        os.path.join(config["results"], "alignment/minimap/from_spoa/spoa.minimap.sam"),

        os.path.join(config["results"], "id_reads/mapped_reads/mapped_seq_id.csv"),
        os.path.join(config["results"], "id_reads/mapped_reads/minimap_output.csv"),
        os.path.join(config["results"], "id_reads/mapped_reads/mapped_consensus.csv"),

        os.path.join(config["results"], "id_reads/filter_id_reads/withinDivergence.csv"),
        os.path.join(config["results"], "id_reads/filter_id_reads/outsideDivergence.csv"),
        os.path.join(config["results"], "id_reads/filter_id_reads/nanDivergence.csv"),

        os.path.join(config["results"], "id_reads/OTU/withinDivergenceOTU.csv"),
        os.path.join(config["results"], "id_reads/OTU/outsideDivergenceOTU.csv"),
        os.path.join(config["results"], "id_reads/OTU/nanDivergenceOTU.csv"),

        os.path.join(config["results"], "id_reads/simple_mapped_reads/simpleMappedWithinDivergence.csv"),
        os.path.join(config["results"], "id_reads/simple_mapped_reads/simpleMappedOutsideDivergence.csv"),
        os.path.join(config["results"], "id_reads/simple_mapped_reads/simpleMappedNaNDivergence.csv"),

        os.path.join(config["results"], "id_reads/cluster_summary/clusterSummaryWithinDivergence.csv"),
        os.path.join(config["results"], "id_reads/cluster_summary/clusterSummaryOutsideDivergence.csv"),
        os.path.join(config["results"], "id_reads/cluster_summary/clusterSummaryNaNDivergence.csv"),

        os.path.join(config["results"], "count_reads/count.reads.barcode.csv"),
        os.path.join(config["results"], "count_reads/count.reads.cutadapt.csv"),
        os.path.join(config["results"], "count_reads/count.reads.filter.csv"),
        os.path.join(config["results"], "count_reads/count.reads.mapping.csv"),

        os.path.join(config["results"], "visuals/plotly/histograms/plotly.barcode.histogram.html"),
        os.path.join(config["results"], "visuals/plotly/histograms/plotly.cutadapt.histogram.html"),
        os.path.join(config["results"], "visuals/plotly/histograms/plotly.filtering.histogram.html"),
        os.path.join(config["results"], "visuals/plotly/histograms/plotly.mapping.histogram.html"),
        os.path.join(config["results"], "visuals/plotly/plotly.box.whisker.html")


"""
Request 2 NVIDIA GPUs, and pass them both to guppy_basecaller
"""
if config["basecall"]["perform_basecall"]:
    checkpoint basecall:
        input: config["basecall_files"]
        output:
            data=directory(os.path.join(config["results"], "basecall")),
            complete=touch(os.path.join(config["results"], ".temp/complete/basecall.complete"))
        params:
            temp_output=os.path.join(config["results"], ".temp/basecall"),
            config=config["basecall"]["configuration"]
        container: config["guppy_container"]
        resources: nvidia_gpu=2
        shell:
            r"""
            command="guppy_basecaller \
            --config {params.config} \
            --input_path {input} \
            --save_path {params.temp_output} \
            --device 'cuda:0,1'"

            # try to resume basecalling. If this does not work, remove the output and try normally
            eval "$command --resume || (rm -rf {params.temp_output} && $command)"


            mv {params.temp_output} {output.data}
            """


    def collate_basecall_input(wildcards):
        checkpoint_output = checkpoints.basecall.get(**wildcards).output
        return glob.glob(os.path.join(checkpoint_output[0], "*.fastq"))
    rule collate_basecall:
        input: collate_basecall_input
        output: fastq_gz=temp(os.path.join(config["results"], ".temp/basecall.merged.files.fastq"))

        shell:
            r"""
            for file in {input}; do
                cat "$file" >> {output.fastq_gz}
            done
            """

    rule NanoQCBasecall:
        input: rules.collate_basecall.output.fastq_gz
        output: directory(os.path.join(config["results"], "visuals/nanoqc/basecall/"))

        shell: "nanoQC -o {output} {input}"

    rule NanoPlotBasecall:
        input: rules.collate_basecall.output.fastq_gz
        output: directory(os.path.join(config["results"], "visuals/nanoplot/basecall/"))

        shell:
            r"""
            NanoPlot --fastq {input} -o {output}
            """


checkpoint barcode:
    input: os.path.join(config["results"], "basecall")
    output: complete=touch(os.path.join(config["results"], ".temp/complete/barcode.complete"))
    params:
        data=temp(directory(os.path.join(config["results"], ".temp/barcode"))),
        guppy_container=config["guppy_container"],
        barcode_kit=config["barcode"]["kit"]
    shell:
        r"""
        guppy_barcoder \
        --input_path {input} \
        --save_path {params.data} \
        --barcode_kits {params.barcode_kit} \
        --recursive
        """

rule merge_barcodes:
    input: lambda wildcards: glob.glob(os.path.join(config["results"], ".temp", "barcode", wildcards.barcode, ".fastq"))
    output: merged=os.path.join(config["results"], "barcode/{barcode}.merged.fastq")
    shell: "cat {input} > {output}"

# cutadapt
rule trim:
    input: rules.merge_barcodes.output.merged
    output: trimmed=os.path.join(config["results"], "trim/{barcode}.trim.fastq")
    params:
        three_prime_adapter=config["cutadapt"]["three_prime_adapter"],
        five_prime_adapter=config["cutadapt"]["five_prime_adapter"],
        error_rate=config["cutadapt"]["error_rate"]
    shell:
        r"""
        cutadapt \
        --revcomp \
        --quiet \
        --adapter {params.three_prime_adapter} \
        --front {params.five_prime_adapter} \
        --error-rate {params.error_rate} \
        --output {output.trimmed} \
        {input}
        """

# NanoFilt
checkpoint filter:
    input: rules.trim.output.trimmed
    output: filter=os.path.join(config["results"], "filter/{barcode}.filter.fastq")
    params:
        min_quality=config["nanofilt"]["min_quality"],
        min_length=config["nanofilt"]["min_filter"],
        max_length=config["nanofilt"]["max_filter"]
    shell:
        r"""
        touch {output}

        NanoFilt \
        --quality {params.min_quality} \
        --length {params.min_length} \
        --maxlength {params.max_length} \
        {input} > {output.filter}
        """

rule filter_gather:
    input: lambda wildcards: expand(os.path.join(config["results"], "filter/{barcode}.filter.fastq"), barcode=return_barcodes(wildcards))
    output: temp(os.path.join(config["results"], ".temp/merge.filter.fastq"))
    shell: "cat {input} > {output}"

def barcode_class_unclass_gather_input(wildcards):
    barcodes = return_barcodes(wildcards)
    return expand(os.path.join(config["results"], "barcode/{barcode}.merged.fastq"),barcode=barcodes)
rule barcode_class_unclass_gather:
    input: barcode_class_unclass_gather_input
    output:
        classified=temp(os.path.join(config["results"], ".temp/barcode.merged.classified.fastq")),
        unclassified=temp(os.path.join(config["results"], ".temp/barcode.merged.unclassified.fastq"))
    shell:
        r"""
        for file in {input}; do
            if [[ "$file" =~ barcode[0-9]{{2}} ]]; then
                cat "$file" >> {output.classified}
            elif [[ "$file" =~ unclassified ]]; then
                cat "$file" >> {output.unclassified}
            fi
        done
        """

rule NanoQCBarcode:
    input:
        classified=rules.barcode_class_unclass_gather.output.classified,
        unclassified=rules.barcode_class_unclass_gather.output.unclassified
    output:
        classified=directory(os.path.join(config["results"], "visuals/nanoqc/barcode/classified")),
        unclassified=directory(os.path.join(config["results"], "visuals/nanoqc/barcode/unclassified"))
    shell:
        r"""
        nanoQC -o {output.classified} {input.classified}
        nanoQC -o {output.unclassified} {input.unclassified}
        """
rule NanoPlotBarcode:
    input:
        classified=rules.barcode_class_unclass_gather.output.classified,
        unclassified=rules.barcode_class_unclass_gather.output.unclassified
    output:
        classified=directory(os.path.join(config["results"], "visuals/nanoplot/barcode/classified")),
        unclassified=directory(os.path.join(config["results"], "visuals/nanoplot/barcode/unclassified"))
    shell:
        r"""
        NanoPlot --fastq {input.classified} -o {output.classified}
        NanoPlot --fastq {input.unclassified} -o {output.unclassified}
        """

rule isONClustPipeline:
    input: rules.filter_gather.output
    output:
        data=directory(os.path.join(config["results"], "isONclust", "pipeline")),
        rule_complete=touch(os.path.join(config["results"], ".temp", "complete", "isONClustPipeline.complete"))
    params:
        aligned_threshold=config["isONclust"]["aligned_threshold"],
        min_fraction=config["isONclust"]["min_fraction"],
        mapped_threshold=config["isONclust"]["mapped_threshold"],
        min_shared = config["isONclust"]["min_shared"]
    shell:
        r"""
        isONclust --ont \
        --fastq {input} \
        --aligned_threshold {params.aligned_threshold} \
        --min_fraction {params.min_fraction} \
        --mapped_threshold {params.mapped_threshold} \
        --min_shared {params.min_shared} \
        --outfolder {output.data}
        """
checkpoint isONclustClusterFastq:
    input:
        pipeline_output=rules.isONClustPipeline.output.data,
        isONClustComplete=rules.isONClustPipeline.output.rule_complete,
        merged_filtering_reads=rules.filter_gather.output
    output:
        cluster_output=directory(os.path.join(config["results"], "isONclust/cluster_fastq")),
        rule_complete=touch(os.path.join(config["results"], ".temp/complete/isONclust.cluster.fastq.complete"))
    params:
        min_quality=config["nanofilt"]["min_quality"],# use same quality as NanoFilt (i.e. rule filtering)
        temp_fastq_input=os.path.join(config["results"], ".temp/merge.filtering.isONclustClusterFastq.fastq")
    shell:
        r"""
        isONclust \
        --q "{params.min_quality}" \
        write_fastq \
        --fastq {input.merged_filtering_reads} \
        --outfolder "{output.cluster_output}" \
        --clusters "{input.pipeline_output}/final_clusters.tsv"
        """


def move_low_reads_input(wildcards):
    """
    We are filtering the *.fastq files from checkpoint.isONclustClusterFastq in this function
    Files that have fewer reads than config["cluster"]["min_reads_per_cluster"] will be added to files_to_move
    These files will be returned in a list to checkpoint move_low_reads
    :param wildcards:
    :return:
    """
    checkpoint_output = checkpoints.isONclustClusterFastq.get(**wildcards).output
    files_to_move = set()
    for file in os.scandir(checkpoint_output[0]):
        if ".fastq" in file.name:

            lines_in_file = open(file.path,"r").readlines()
            count_lines_in_file = len(lines_in_file)
            count_reads_in_file = count_lines_in_file / 4

            # keep files that are equal to OR GREATER than the cutoff
            if count_reads_in_file < config["cluster"]["min_reads_per_cluster"]:
                # only get the file name (remove the extension)
                files_to_move.add(file.name.split(".")[0])
    return expand(os.path.join(checkpoint_output[0], "{file_move}.fastq"), file_move=files_to_move)
checkpoint move_low_reads:
    input: move_low_reads_input
    output:
        data=directory(os.path.join(config["results"], "LowClusterReads")),
        complete=touch(os.path.join(config["results"], ".temp/complete/remove.low.reads.complete"))
    shell:
        r"""
        mkdir -p {output.data}
        for file_path in {input}; do
            file_name=$(basename -- "$file_path")
            mv "$file_path" {output.data}/"$file_name"
        done
        """


def spoa_input(wildcards):
    isonclust_output = checkpoints.isONclustClusterFastq.get(**wildcards).output
    move_low_reads_output = checkpoints.move_low_reads.get(**wildcards).output
    return glob.glob(os.path.join(isonclust_output[0], "*.fastq"))
rule spoa:
    input: spoa_input
    output: os.path.join(config["results"], "spoa/consensus.sequences.fasta")
    params: temp_spoa=os.path.join(config["results"], ".temp/spoa.temp.fasta")
    run:
        os.makedirs(os.path.join(config["results"], "spoa"), exist_ok=True)

        # remove temp output (in case it exists from previous run)
        shell("rm -f {params.temp_spoa}")

        for file in str(input).split(" "):

            # perform spoa
            shell("spoa {input} -r 0 > {params.temp_spoa}")

            # get file name without extension (i.e. /path/to/file/35.fasta -> 35)
            basename = os.path.basename(file)
            cluster_number = basename.split(".")[0]

            # overwrite first line in file with `>cluster_{cluster_number}`
            file_lines = open(str(params.temp_spoa),"r").readlines()
            file_lines[0] = f">cluster_{cluster_number}\n"

            # append new lines to output file
            open(str(output),"a").writelines(file_lines)

        # remove temp output, it is no longer needed
        shell("rm -f {params.temp_spoa}")


def minimap_from_filtering_input(wildcards):
    checkpoint_output = checkpoints.filter.get(**wildcards).output
    return glob.glob(os.path.join(config["results"], f"filter/{wildcards.barcode}.filter.fastq"))
rule minimap_from_filtering:
    input: minimap_from_filtering_input
    output: os.path.join(config["results"], "alignment/minimap/from_filtering/{barcode}.minimap.sam")
    params: alignment_reference=config["reference_database"]
    shell:
        r"""
        minimap2 \
        -ax map-ont \
        {params.alignment_reference} \
        {input} > {output}
        """

rule minimap_from_spoa:
    input: rules.spoa.output[0]
    output: os.path.join(config["results"], "alignment/minimap/from_spoa/spoa.minimap.sam")
    params: alignment_reference=config["reference_database"]
    shell:
        r"""
        minimap2 \
        -ax map-ont \
        {params.alignment_reference} \
        {input} > {output}
        """


# TODO: It appears that mapped_consensus_csv is missing a header for the first column
rule id_reads:
    input:
        filtering=lambda wildcards: expand(os.path.join(config["results"], "filter", "{barcode}.filter.fastq"), barcode=return_barcodes(wildcards)),
        clustering=rules.isONClustPipeline.output[0],
        minimap=rules.minimap_from_spoa.output[0]
        # filtering=filtering_output
    output:
        mapped_seq_id_csv=os.path.join(config["results"], "id_reads/mapped_reads/mapped_seq_id.csv"),
        minimap_output=os.path.join(config["results"], "id_reads/mapped_reads/minimap_output.csv"),
        mapped_consensus_csv=os.path.join(config["results"], "id_reads/mapped_reads/mapped_consensus.csv")
    params: results_folder=config["results"]
    script: "scripts/id_reads.py"

rule filter_id_reads_mapped_sequence:
    input: csv=rules.id_reads.output.mapped_seq_id_csv
    output:
        within_divergence=os.path.join(config["results"], "id_reads/filter_id_reads/withinDivergence.csv"),
        outside_divergence=os.path.join(config["results"], "id_reads/filter_id_reads/outsideDivergence.csv"),
        nan_divergence=os.path.join(config["results"], "id_reads/filter_id_reads/nanDivergence.csv")
    params: divergence_threshold=config["cluster"]["divergence_threshold"]
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
        within_divergence_otu=os.path.join(config["results"], "id_reads/OTU/withinDivergenceOTU.csv"),
        outside_divergence_otu=os.path.join(config["results"], "id_reads/OTU/outsideDivergenceOTU.csv"),
        nan_divergence_otu=os.path.join(config["results"], "id_reads/OTU/nanDivergenceOTU.csv")
    script: "scripts/generateOTU.py"

rule simple_mapped_sequence_id:
    input: rules.id_reads.output.mapped_seq_id_csv
    output:
        within_divergence=os.path.join(config["results"], "id_reads/simple_mapped_reads/simpleMappedWithinDivergence.csv"),
        outside_divergence=os.path.join(config["results"], "id_reads/simple_mapped_reads/simpleMappedOutsideDivergence.csv"),
        nan_divergence=os.path.join(config["results"], "id_reads/simple_mapped_reads/simpleMappedNaNDivergence.csv")
    params: divergence_threshold=config["cluster"]["divergence_threshold"]
    script: "scripts/simpleMappedSequenceID.py"

rule cluster_summary:
    input: rules.id_reads.output.mapped_seq_id_csv
    output:
        within_divergence=os.path.join(config["results"], "id_reads/cluster_summary/clusterSummaryWithinDivergence.csv"),
        outside_divergence=os.path.join(config["results"], "id_reads/cluster_summary/clusterSummaryOutsideDivergence.csv"),
        nan_divergence=os.path.join(config["results"], "id_reads/cluster_summary/clusterSummaryNaNDivergence.csv")
    params: divergence_threshold=config["cluster"]["divergence_threshold"]
    script: "scripts/clusterSummary.py"


rule count_reads_barcode:
    input:
        lambda wildcards: expand(os.path.join(config["results"], "barcode", "{barcode}.merged.fastq"), barcode=return_barcodes(wildcards))
        #count_reads_barcode_input
    output: os.path.join(config["results"], "count_reads/count.reads.barcode.csv")
    params: process="barcode"
    script: "scripts/CountReads.py"


rule count_reads_cutadapt:
    input:
        #count_reads_cutadapt_input
        lambda wildcards: expand(os.path.join(config["results"], "trim/{barcode}.trim.fastq"), barcode=return_barcodes(wildcards))
    output: os.path.join(config["results"], "count_reads/count.reads.cutadapt.csv")
    params: process="cutadapt"
    script: "scripts/CountReads.py"


rule count_filtering:
    input:
        #count_filtering_input
        lambda wildcards: expand(os.path.join(config["results"], "filter/{barcode}.filter.fastq"), barcode=return_barcodes(wildcards))
    output: os.path.join(config["results"], "count_reads/count.reads.filter.csv")
    params: process="filtering"
    script: "scripts/CountReads.py"


rule count_reads_mapping:
    input:
        #count_minimap_reads
        lambda wildcards: expand(os.path.join(config["results"], "alignment/minimap/from_filtering/{barcode}.minimap.sam"),
            barcode=return_barcodes(wildcards))
    output:os.path.join(config["results"], "count_reads/count.reads.mapping.csv")
    params: process="mapping"
    script: "scripts/CountReads.py"


rule plotly_barcode_histogram:
    input: rules.count_reads_barcode.output[0]
    output: os.path.join(config["results"], "visuals/plotly/histograms/plotly.barcode.histogram.html")
    params: sub_title="Performed after Merging Files"
    script: "scripts/PlotlyHistogram.py"

rule plotly_cutadapt_histogram:
    input: rules.count_reads_cutadapt.output[0]
    output: os.path.join(config["results"], "visuals/plotly/histograms/plotly.cutadapt.histogram.html")
    params: sub_title="Performed after Cutadapt"
    script: "scripts/PlotlyHistogram.py"

rule plotly_filtering_histogram:
    input: rules.count_filtering.output[0]
    output: os.path.join(config["results"], "visuals/plotly/histograms/plotly.filtering.histogram.html")
    params: sub_title="Performed after Filtering"
    script: "scripts/PlotlyHistogram.py"

rule plotly_mapping_histogram:
    input: rules.count_reads_mapping.output[0]
    output: os.path.join(config["results"], "visuals/plotly/histograms/plotly.mapping.histogram.html")
    params: sub_title="Performed after Mapping"
    script: "scripts/PlotlyHistogram.py"

rule plotly_box_whisker_generation:
    input:
        rules.count_reads_barcode.output[0],
        rules.count_reads_cutadapt.output[0],
        rules.count_filtering.output[0],
        rules.count_reads_mapping.output[0]
    output: os.path.join(config["results"], "visuals/plotly/plotly.box.whisker.html")
    script: "scripts/PlotlyBoxWhisker.py"
