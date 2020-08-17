import os
from pathlib import Path
import glob
configfile: "config.yaml"

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
    basecall_checkpoint = checkpoints.basecall.get(**wildcards).output[0]
    barcode_checkpoint = checkpoints.barcode.get(**wildcards).output[0]
    barcodes = set()  # a set is like a list, but only stores unique values
    for folder in os.listdir(barcode_checkpoint):
        full_path = os.path.join(barcode_checkpoint, folder)
        if Path(full_path).is_dir():
            barcodes.add(folder)

    merge_files = [config['results_folder'] + "barcode/" + barcode + ".merged.fastq" for barcode in barcodes]
    return merge_files
def nanoqc_basecall_data(wildcards):
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return config['results_folder'] + "visuals/nanoqc/basecall/"
def nanoqc_barcode_classified(wildcards):
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return config['results_folder'] + "visuals/nanoqc/barcode/classified"
def nanoqc_barcode_unclassified(wildcards):
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return config['results_folder'] + "visuals/nanoqc/barcode/unclassified"
def cutadapt(wildcards):
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return expand(
        config['results_folder'] + "cutadapt/{barcode}.cutadapt.fastq",
        barcode=return_barcode_numbers(checkpoint_output))
def filtering(wildcards):
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return expand(
        config['results_folder'] + "filter/{barcode}.filter.fastq",
        barcode=return_barcode_numbers(checkpoint_output))
def guppy_aligner(wildcards):
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return expand(
        config['results_folder'] + "alignment/guppy/sam_files/{barcode}.guppy.sam",
               barcode=return_barcode_numbers(checkpoint_output))
def minimap_aligner(wildcards):
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return expand(
        config['results_folder'] + "alignment/minimap/{barcode}.minimap.sam",
        barcode=return_barcode_numbers(checkpoint_output))
def vsearch_aligner(wildcards):
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return expand(
        config['results_folder'] + "alignment/vsearch/{barcode}.vsearch.uc",
        barcode=return_barcode_numbers(checkpoint_output))
def nanoplot_basecall(wildcards):
    checkpoint_output = checkpoints.basecall.get(**wildcards).output[0]
    return config['results_folder'] + "visuals/nanoplot/basecall/"
def nanoplot_barcode_classified(wildcards):
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return config['results_folder'] + "visuals/nanoplot/barcode/classified"
def nanoplot_barcode_unclassified(wildcards):
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return config['results_folder'] + "visuals/nanoplot/barcode/unclassified"
def plotly_histogram_barcode(wildcards):
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return config['results_folder'] + "visuals/plotly/histograms/plotly.barcode.histogram.html"
def plotly_histogram_cutadapt(wildcards):
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return config['results_folder'] + "visuals/plotly/histograms/plotly.cutadapt.histogram.html"
def plotly_histogram_filtering(wildcards):
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return config['results_folder'] + "visuals/plotly/histograms/plotly.filtering.histogram.html"
def plotly_histogram_mapping(wildcards):
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return config['results_folder'] + "visuals/plotly/histograms/plotly.mapping.histogram.html"
def plotly_box_whisker(wildcards):
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return config['results_folder'] + "visuals/plotly/plotly.box.whisker.html"
FAST5_FILES = glob_wildcards(config['fast5_location'] + "{fast5_file}.fast5").fast5_file



rule all:
    input:
        barcode_merge_files,  #......................................... Basecall, barcode, and merge files in a checkpoint
        nanoqc_basecall_data,  #....................................... NanoQC after basecall
        nanoqc_barcode_classified,  #.................................. NanoQC classified barcodes
        nanoqc_barcode_unclassified,  #................................ NanoQC unclassified barcodes
        cutadapt,  #................................................... Trim (Cutadapt)
        filtering,  #................................................... NanoFilt
        guppy_aligner,  #.............................................. Guppy Aligner
        minimap_aligner,  #............................................ MiniMap Aligner
        vsearch_aligner,  #............................................ VSearch Aligner
        nanoplot_basecall,  #................................................... NanoPlot
        nanoplot_barcode_classified,
        nanoplot_barcode_unclassified,
        plotly_histogram_barcode,  #................................... Plotly barcode histogram
        plotly_histogram_cutadapt,  #.................................. Plotly cutadapt histogram
        plotly_histogram_filtering,  #................................. Plotly filtering histogram
        plotly_histogram_mapping,  #................................... Plotly mapping histogram
        plotly_box_whisker  #.......................................... Plotly box and whisker plot



checkpoint basecall:
    input:
        config['fast5_location']
    output:
        output = directory(config['results_folder'] + "basecall/")
    params:
        configuration = config["basecall_configuration"],
        callers = config['num_callers'],
        threads_per_caller = config['num_threads_per_caller']
    shell:
        r"""
        echo Basecalling
        
        guppy_basecaller \
        --config {params.configuration} \
        --input_path {input} \
        --save_path {output} \
        --num_callers {params.callers} \
        --cpu_threads_per_caller {params.threads_per_caller} \
        --recursive
        """



checkpoint barcode:
    input:
        rules.basecall.output[0]
    output:
        directory(config['results_folder'] + ".temp/barcodeTempOutput/")
    params:
        barcode_kit = config['barcode_kit']
    shell:
        r"""
        echo Barcoding
        
        guppy_barcoder \
        -i {input} \
        -s {output} \
        --barcode_kits {params.barcode_kit} \
        --recursive \
        --quiet       
        """



def merge_files_input(wildcards):
    import glob
    return glob.glob(config['results_folder'] + f".temp/barcodeTempOutput/{wildcards.barcode}/*.fastq")
rule merge_files:
    input:
        merge_files_input
    output:
        config['results_folder'] + "barcode/{barcode}.merged.fastq"
    params:
        input_folder = config['results_folder'] + ".barcodeTempOutput",
        save_folder = config['results_folder'] + "barcode"
    shell:
        r"""
        for item in {input}; do
            cat $item >> {output}
        done
        """



def create_classified_unclassified_input(wildcards):
    checkpoint_output = checkpoints.barcode.get(**wildcards).output[0]
    return glob.glob(config['results_folder'] + ".temp/barcodeTempOutput/**/*.fastq")
rule collate_basecall_fastq_files:
    input:
        config['results_folder'] + "basecall"
    output:
        config['results_folder'] + ".temp/basecall.temp.merged.files.fastq"
    shell:
        r"""        
        # concatenate each file in the params directory to the output file
        for file in "{input[0]}"/*.fastq; do
            cat "$file" >> {output}
        done
        """
rule create_classified_unclassified_file:
    input:
        create_classified_unclassified_input
    output:
        classified = config['results_folder'] + ".temp/barcode.classified.merged.temp.fastq",
        unclassified = config['results_folder'] + ".temp/barcode.unclassified.merged.temp.fastq"
    shell:
        r"""
        touch {output.classified}
        touch {output.unclassified}
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
        directory(config['results_folder'] + "visuals/nanoqc/basecall/")
    shell:
        r"""
        nanoQC -o {output} {input}
        """
rule NanoQCBarcode:
    input:
        classified = rules.create_classified_unclassified_file.output.classified,
        unclassified = rules.create_classified_unclassified_file.output.unclassified
    output:
        classified = directory(config['results_folder'] + "visuals/nanoqc/barcode/classified"),
        unclassified = directory(config['results_folder'] + "visuals/nanoqc/barcode/unclassified")
    shell:
        r"""
        nanoQC -o {output.classified} {input.classified}
        nanoQC -o {output.unclassified} {input.unclassified}
        """



rule NanoPlotBasecall:
    input:
        rules.collate_basecall_fastq_files.output[0]
    output:
        directory(config['results_folder'] + "visuals/nanoplot/basecall/")
    shell:
        r"""
        NanoPlot --fastq {input} -o {output}
        """
rule NanoPlotBarcode:
    input:
        classified = rules.create_classified_unclassified_file.output.classified,
        unclassified = rules.create_classified_unclassified_file.output.unclassified

    output:
        classified = directory(config['results_folder'] + "visuals/nanoplot/barcode/classified"),
        unclassified = directory(config['results_folder'] + "visuals/nanoplot/barcode/unclassified")
    shell:
        r"""
        NanoPlot --fastq {input.classified} -o {output.classified}
        NanoPlot --fastq {input.unclassified} -o {output.unclassified}
        """


rule cutadapt:
    input:
        rules.merge_files.output[0]
    output:
        config['results_folder'] + "cutadapt/{barcode}.cutadapt.fastq"
    params:
        three_prime_adapter = config['trim_three_prime_adapter'],
        five_prime_adapter = config['trim_five_prime_adapter'],
        error_rate = config['trim_error_rate']
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



rule filtering:
    input:
        rules.cutadapt.output[0]
    output:
        config['results_folder'] + "filter/{barcode}.filter.fastq"
    shell:
        r"""
        touch {output}
        NanoFilt {input} > {output}
        """



rule guppy_aligner:
    input:
        rules.filtering.output[0]
    output:
        sam_files = config['results_folder'] + "alignment/guppy/sam_files/{barcode}.guppy.sam",
        alignment_summary = config['results_folder'] + "alignment/guppy/alignment_summary/{barcode}.alignment.summary.csv",
        log_file = config['results_folder'] + "alignment/guppy/logs/{barcode}.guppy.log",
    params:
        barcode = "{barcode}",
        temp_dir = config['results_folder'] + ".temp/guppy",
        alignment_reference = config['alignment_reference_file']
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
rule minimap_aligner:
    input:
        rules.filtering.output[0]
    output:
        config['results_folder'] + "alignment/minimap/{barcode}.minimap.sam"
    params:
        alignment_reference = config['alignment_reference_file']
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
        rules.filtering.output[0]
    output:
        temp(config['results_folder'] + "alignment/vsearch/.{barcode}.temp.fasta")
    shell:
        r"""
        seqkit fq2fa {input} > {output}
        """
rule vsearch_aligner:
    input:
        rules.fq2fa.output[0]
    output:
        config['results_folder'] + "alignment/vsearch/{barcode}.vsearch.uc"
    params:
        alignment_reference = config['alignment_reference_file']
    shell:
        r"""
        vsearch \
        --usearch_global \
        {input} \
        --db {params.alignment_reference} \
        --id 0 \
        --uc {output} \
        --quiet
        """




rule count_reads_barcode:
    input:
        expand(rules.merge_files.output[0],
               barcode=glob_wildcards(config['results_folder'] + "barcode/{barcode}.merged.fastq").barcode)
    output:
        config['results_folder'] + "count_reads/count.reads.barcode.csv"
    params:
        process = "barcode"
    script:
        "scripts/CountReads.py"
rule count_reads_cutadapt:
    input:
        expand(rules.cutadapt.output[0],
               barcode=glob_wildcards(config['results_folder'] + "cutadapt/{barcode}.cutadapt.fastq").barcode)
    output:
        config['results_folder'] + "count_reads/count.reads.cutadapt.csv"
    params:
        process = "cutadapt"
    script:
        "scripts/CountReads.py"
rule count_reads_filtering:
    input:
        expand(rules.filtering.output[0],
               barcode=glob_wildcards(config['results_folder'] + "filter/{barcode}.filter.fastq").barcode)
    output:
        config['results_folder'] + "count_reads/count.reads.filter.csv"
    params:
        process = "filtering"
    script:
        "scripts/CountReads.py"
rule count_reads_mapping:
    input:
        expand(rules.guppy_aligner.output.alignment_summary,
               barcode=glob_wildcards(config['results_folder'] + "alignment/guppy/alignment_summary/{barcode}.alignment.summary.csv").barcode)
    output:
        config['results_folder'] + "count_reads/count.reads.mapping.csv"
    params:
        process = "mapping"
    script:
        "scripts/CountReads.py"



rule plotly_barcode_histogram:
    input:
        rules.count_reads_barcode.output[0]
    output:
        config['results_folder'] + "visuals/plotly/histograms/plotly.barcode.histogram.html"
    params:
        sub_title = "Performed after Merging Files"
    script:
        "scripts/PlotlyHistogram.py"
rule plotly_cutadapt_histogram:
    input:
        rules.count_reads_cutadapt.output[0]
    output:
        config['results_folder'] + "visuals/plotly/histograms/plotly.cutadapt.histogram.html"
    params:
        sub_title = "Performed after Cutadapt",
    script:
        "scripts/PlotlyHistogram.py"
rule plotly_filtering_histogram:
    input:
        rules.count_reads_filtering.output[0]
    output:
        config['results_folder'] + "visuals/plotly/histograms/plotly.filtering.histogram.html"
    params:
        sub_title = "Performed after Filtering",
    script:
        "scripts/PlotlyHistogram.py"
rule plotly_mapping_histogram:
    input:
        rules.count_reads_mapping.output[0]
    output:
        config['results_folder'] + "visuals/plotly/histograms/plotly.mapping.histogram.html"
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
        config['results_folder'] + "visuals/plotly/plotly.box.whisker.html"
    script:
        "scripts/PlotlyBoxWhisker.py"

