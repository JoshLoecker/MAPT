# Pipeline (MAPT)

## Naming
I would like to name this project as it has been my first, and most in-depth, project to date. 
Due to the nature of this project revolving around mapping noisy microbial reads (hopefully) down to the species level, 
the following name seemed appropriate: Microbial Automated Processing Tool, or MAPT for short.

### Installation and Running
This project is meant to run using Singularity or Docker
### Folder Setup
Two paths are needed with this container: 1) A `results` folder, and 2) A `data` folder. These folders can be named as you please, but this guide will use these respective names.  

Under the `data` folder, a folder **must** be named `fast5` **or** `fastq`. If this is not done, guppy_basecaller/barcoder will not be able to find any 
`.fast5` files, even if they are present, as it looks for `.fast5` files under the `fast5` folder.
The `results` folder must exist before running the container. Your alignment file should also be placed within the `data` folder.
Your `data` folder structure may look as follows:
```
home
| -- Rob
    | -- Projects
        | -- cache_dir
        | -- alignment_files
            | -- silva_alignment_file.fasta
            | -- another_alignment_file.fasta
        | -- run_1
            | -- data
                | -- fast5
                    | -- file_1.fast5
                    | -- file_2.fast5
                    | -- file_3.fast5
                | -- fastq
                    | -- file_1.fastq
                    | -- file_2.fastq
                    | -- file_3.fastq
                | -- some_other_file.txt
            | -- results
```

### To Start
Set up a few variables in the terminal, using the above structure as an example. 
Change these values for your workflow.

    results="/home/Rob/Projects/run_1/results"
    data="/home/Rob/Projects/run_1/data/"
    alignment_path="/home/Rob/Projects/alignment_files/silva_alignment_file.fasta"
    num_basecallers=3
    threads_per_caller=5
    basecall=True  # True or False, should basecalling be done?
    basecall_configuration=""
    barcode_kit=""
    cutadapt_trim_error_rate=0.15  # default value of 0.15
    cutadapt_trim_three_prime_adapter=""
    cutadapt_trim_five_prime_adapter=""
    cluster_cutoff=3  # Clusters with this many reads OR GREATER are kept. Clusters with fewer are moved to /results/.temp/TooFewReadsInCluster
    mapped_reads_divergence_threshold=0.05  # The divergence threshold for mapped reads after filtering & clustering
    nanofilt_filtering_min=1000  # default value of 1000
    nanofilt_filtering_max=2000  # default value of 2000
    

   If you do not have fast5 files (i.e. no need for basecalling), set the `perform_basecall` option to `False`.
  
   
### Singularty
1. If you chose to work with Singularity, ensure it is already installed on your system by running `singularity --version`. 
   If `singularity version . . .` does not appear, [install singularity here](https://singularity.lbl.gov/install-linux). 
   Singularity is already installed on most clusters, such as SciNET

2. Next, pull the docker container  
    `singularity pull docker://joshloecker/pipeline:latest`

4. Assuming the `To Start` section was completed, the following command
can safely be copied and pasted. This will download the container, but it is
not yet running
    ```
    singularity run \
    --bind $results:/results \
    --bind $data:/data \
    --bind $cache_dir /home/$USER/.cache \
    --env alignment_path=$alignment_path \
    --env num_basecallers=$num_basecallers \
    --env threads_per_caller=$threads_per_caller \
    --env basecall_configuration=$basecall_configuration \
    --env barcode_kit=$barcode_kit \
    --env cutadapt_trim_error_rate=$cutadapt_trim_error \
    --env cutadapt_trim_three_prime_adapter=$cutadapt_trim_three_prime_adapter \
    --env cutadapt_trim_five_prime_adapter=$cutadapt_trim_five_prime_adapter \
    --env nanofilt_filtering_min=$nanofilt_filtering_min \
    --env nanofilt_filtering_max=$nanofilt_filtering_max \
    --env cluster_cutoff=$cluster_cutoff \
    --env mapped_reads_divergence_threshold=$mapped_reads_divergence_threshold \
    --env basecall=$basecall \
    docker://joshloecker/pipeline:latest
	```
   
    This will run the container

4. If you would like to see a dry-run of the pipeline, append `--dry-run`
to the end of the `run` command. [Other Snakemake flags](https://snakemake.readthedocs.io/en/stable/executing/cli.html)
can also be added in this manner.  
   ```
   singularity run joshloecker/pipeline --dry-run
   ```
    Note: Due to this workflow running in a container, not all flags have been tested
    nor are confirmed to work as expected.  
   
   
### Docker
1. If you chose to work with Docker, ensure Docker is already installed on your 
system by running `docker --version`. If `Docker version . . .` does not appear,
[install docker here](https://docs.docker.com/get-docker/)  

2. If you would like to start by just downloading the container, run the following  
    ```docker pull joshloecker/pipeline:latest```

2. Then, run the container using the following command. This can safely be copied
and pasted, assuming the `To Start` step has been completed
    ```
    docker run -it \
    --name=pipeline \
    --mount type=bind,source=$results,target=/results \
    --mount type=bind,source=$data,target=/data \
    joshloecker/pipeline \
    --config \
    alignment_path=$alignment_path \
    num_basecallers=$num_basecallers \
    threads_per_caller=$threads_per_caller \
    basecall_configuration=$basecall_configuration \
    barcode_kit=$barcode_kit \
    cutadapt_trim_error_rate=$cutadapt_trim_error_rate \
    cutadapt_trim_three_prime_adapter=$cutadapt_trim_three_prime_adapter \
    cutadapt_trim_five_prime_adapter=$cutadapt_trim_five_prime_adapter \
    nanofilt_filtering_min=$nanofilt_filtering_min \
    nanofilt_filtering_max=$nanofilt_filtering_max \
    --env cluster_cutoff=$cluster_cutoff \
    --env mapped_reads_divergence_threshold=$mapped_reads_divergence_threshold \
    basecall=$basecall
    ```
    This will start the container, or download it if it is not downloaded.
    <br><br>
    If you would like to do a dry-run before starting the container, simply add `--dry-run` to the end of the command.
    You must delete the container (`docker rm pipeline`) and re-run it without `--dry-run` to begin your workflow.
    Any additional snakemake flags should be able to be entered at the end of 
    this command as well.

### Workflow
The following workflow will be completed, relatively in this order
1. [Oxford Nanopore Basecalling](https://community.nanoporetech.com/protocols/Guppy-protocol/v/gpb_2003_v1_revt_14dec2018) (i.e. guppy_basecaller)
2. [Oxford Nanopore Barcoding](https://community.nanoporetech.com/protocols/Guppy-protocol/v/gpb_2003_v1_revt_14dec2018) (i.e. guppy_barcoder)
3. [NanoQC](https://github.com/wdecoster/nanoQC)  
    a. NanoQC will be called independently on basecalling and barcoding results
4. [NanoPlot](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwiz89ql2cHsAhUjAp0JHUVoCFAQFjAAegQIARAC&url=https%3A%2F%2Fgithub.com%2Fwdecoster%2FNanoPlot&usg=AOvVaw00LEGNovoQzjS5KCUxwD0v)
5. [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)
6. [NanoFilt](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwigp9S92cHsAhWSLc0KHYp8C-oQFjAAegQIBBAC&url=https%3A%2F%2Fgithub.com%2Fwdecoster%2Fnanofilt&usg=AOvVaw20npdGb-VRvmFH1SY6-P6C)
7. [isONclust](https://github.com/ksahlin/isONclust)
8. [SPOA](https://github.com/rvaser/spoa)
9. [Minimap Aligner](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwjUtqXD2cHsAhVCGM0KHYsKAqcQFjAAegQIARAC&url=https%3A%2F%2Fgithub.com%2Flh3%2Fminimap2&usg=AOvVaw3UvK2vgr0fj_4GS68K8V26) on SPOA and Filtering results
10. [IsoCon](https://github.com/ksahlin/IsoCon)
11. [Plotly](https://pypi.org/project/plotly/) Visuals

### Results
Results from the pipeline will be saved to the path you bound during the **Installation and Running** step (i.e. `/home/Rob/Projects/run_1/results`). 
From here, you will be able to see each of the folders that are output, such as `basecalling`, `barcoding`, `visuals`, etc.


### Known Errors
Deleting files from the `.temp` folder will cause the pipeline to regenerate these files, along with any output downstream.  
If you experience any other errors, contact [joshua.loecker@jacks.sdstate.edu](mailto:joshua.loecker@jacks.sdstate.edu) for assistance

<!--stackedit_data:
eyJoaXN0b3J5IjpbOTIyMDg4NzE1XX0=
-->
