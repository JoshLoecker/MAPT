# MAPT

### Naming
This has been my most notable, and in-depth, project to date. As a result, I feel it should have a proper name.
Due to the nature of this project revolving around mapping noisy microbial reads (hopefully) down to the species level, 
the following name seemed appropriate: Microbial Automated Processing Tool, or MAPT for short.


### Usage

This is a workflow/pipeline aimed to streamline post-processing of noisy reads from [Oxford Nanopore MinION Sequencing Data](https://www.nanoporetech.com). The following software will need to be installed in a conda environment. A startup-script is provided in the GitHub repository. This is explained further on this page.

The following software will be used (not necessarily in this order, and this list may be incomplete):
- Guppy Basecalling
- Guppy Barcoding
- Cutadapt
- NanoFilt
- isONclust
- spoa
- NanoPlot
- Plotly Graphs


### Installation & Setup

This project was built and testing on the following versions of Singularity and Docker. Reproducibility and stability cannot be guaranteed on earlier versions.
1. Singularity (version 3.6.4-1.el7 or higher)
2. Docker (version 20.10.0 or higher)

[Download singularity](https://singularity.lbl.gov/)  
[Download docker](https://www.docker.com/products/docker-desktop)

Git and Miniconda are required to download the pipeline and create a new Conda environment with the required software

[Download git](https://git-scm.com/downloads)
[Download miniconda](https://docs.conda.io/en/latest/miniconda.html)

1. Start by cloning the GitHub repo  
    a. `git clone https://github.com/JoshLoecker/MAPT`  
    b. This will create a new folder `MAPT` in your current directory. You may move or rename this folder if you like
    c. While a current version of the pipeline exists at `/project/brookings_minion/pipeline`, it is not possible to have multiple runs starting from the same working directory
    1) Attempting to do so will result in an error, such as `IncompleteFilesException, the files below seem to be incomplete`
	2) This is because snakemake attempts to determine what outputs need to be created, and 'sees' that output is unfinished (from a currently running job, using this directry)
	3) The easiest fix is simply cloning a new version of the pipeline. This can even be done in your home directory, as the pipeline is just 4M in size
   
2. Installing the environment
	a. **Note**: If you are on SciNet, a conda environment already exists. It is located at `/project/brookings_minoin/conda-envs/mapt_pipeline`
	1) Activating the environment can be seen in next step
    
	b. On SciNet, storage space in your home directory is limited. It is [strongly recommended](https://scinet.usda.gov/guide/ceres/#quotas-on-home-and-project-directories) to store environments in your `/project` directory.  
    c. By saving environments in your `/project` directory, it is possible for everyone on your team to use the same environment, even if multiple people are running jobs at the same time.  
    d. To create a new conda environment, execute the following: `conda env create -f environment.yaml -p /path/to/my/data/env_name`  
    1) This will create a new environment with the name `env_name`.
	
3. Activating the environment
   a. To activate the environment, its location must be known.
   1) This assumes there is currently a `mapt_pipeline` environment installed at `/project/brookings_minion/conda-envs/mapt_pipeline`
	
	b. Simply type `conda activate /project/brookings_minion/conda-envs/mapt_pipeline` to activate the environment
    c. By doing so, you now have access to all packages installed within
    d. If you see the full path before your username (i.e., `(/project/brookings_minion/conda-envs/mapt_pipeline) USERNAME /current/directory $ |` please do the following  
	1) `nano ~/.condarc`
	2) Add the following: `env_prompt: ({name})`
	3) Save with `CTRL + X` -> `y` -> `[ENTER]`
	
	e. To deactivate the environment, type `conda deactivate`.
  
4. To be able to use the Guppy suite of tools, a singularity container is available
   a. This container will be used by snakemake during its run.  
   b. For more information, view Running Guppy in a Singularity Container at the bottom of this page
5. The final step is to edit several lines within the `pipeline/config.yaml` file  
    a. First, set the `results`, `basecall_files` (or `barcode_files`), `reference_database`, and `guppy_container` to their appropriate locations  
	1) `results` is where you would like results of the pipeline to be stored.  
	2) `basecall_files` holds your fast5/fastq files for Guppy.
	   a) If you are going to skip basecalling and do only barcoding, the `basecall_files` config can be left empty (i.e. `""`)
	3) `reference_database` is the database you will be using with MiniMap for alignments.  
	4) `guppy_container` is the location of the guppy container you will be using.  
		a) It can most likely be left as-is
	
	b. Set any other values required under the `DEFAULT VALUES` section.

**NOTE**: Any values not changed during in the `pipeline/config.yaml` file will be used during the execution of the pipeline. **This can lead to unexpected outcomes**.

### Running the Pipeline
Once these steps are done, the pipeline is ready to run. The pipeline can be run in several methods

1. Interactive Runs  
    a. If you would like to see the output of jobs as they happen, or you have a short job you would like to ensure is working, Interactive Runs can be useful  
    b. Follow the [guide](https://scinet.usda.gov/guide/ceres/#interactive-mode) for help on how to set up an interactive run
    1) A list of available partitions & queues can be [found here](https://scinet.usda.gov/guide/ceres/#partitions-or-queues)
    1) In short, the following structure should be used: `srun --pty -p [QUEUE_CHOICE] -t hh:mm:ss -n [TASKS] -N [NODES] /bin/bash -l`  
   
	c. It should be noted that if you have an Interactive Run, and your connection to the server is lost, the job will quit immediately. Because of this, it is recommended to use Slurm Jobs instead  
2. Slurm Jobs  
	a. If you would like to close your connection to the server, Slurm Jobs are the most versatile tool.  
    b. [SciNet User Guide](https://scinet.usda.gov/guide/ceres/)  
    1) [Requesting the proper nodes and cores](https://scinet.usda.gov/guide/ceres/#requesting-the-proper-number-of-nodes-and-cores)  
    2) [Simple how-to on SLURM](https://scinet.usda.gov/guide/ceres/#batch-mode)  
	
	c. For more SBATCH options, see [this guide](https://osirim.irit.fr/site/en/articles/sbatch-options)    
	d. For some simple getting-started scripts, view the example slurm scripts in one of the two locations:
    1) SciNet: `/project/brookings_minion/example_slurm_scripts/`
	2) GitHub: [Example SLURM Scripts](https://github.com/JoshLoecker/pipeline/tree/master/Example%20SLURM%20Scripts)

3. Dry Runs  
    a. Dry-Runs can be done in Interactive Runs or Slurm Jobs.  
	b. A dry-run allows you to see what steps need to be done, and ensures preliminary configuration is set up correctly, without actually executing any portion of the workflow.  
    c. To perform a dry-run, activate the `mapt_pipeline` conda environment by running the following command: `conda activate /project/brookings_minion/conda-envs/mapt_pipeline`  
    d. Next, call snakemake with a dry-run: `snakemake -j 1 -n`  
    1) `snakemake`: Start snakemake  
	2) `-j 1` (or `--cores 1`): Use 1 core for the dry-run (we do not want to bog-down the login node on the cluster. More cores will not speed this process up)
	4) `-n` (or `--dry-run`): This is the dry-run flag for snakemake  
	
	d. This will output a fair amount of information, showing what rules need to be completed. The pipeline can then be ran following point 1 or 2 above.  

## Notes to Future Maintainers
### How to build a new Guppy container

Building a new guppy container is relatively simple
1. Go to [SyLabs Cloud Builder](https://cloud.sylabs.io/builder)
	a. Sign in with one of the options available
	b. Click on your username in the top right corner
	c. On the drop down, click "Access Tokens"
2. Generate a new Access Token
	a. Copy this token, and save it in a secure location. It will not be able to be accessed again
3. [Log in to SciNet](https://scinet.usda.gov/guide/ceres/#system-access)
4. Navigate to /project/brookings_minion/ on SciNet
	a. There is most likely a file named `guppy_container.sif`
	b. We are going to be updating this file
5. Run `singularity remote login` and follow the instructions to enable the remote build server
6. Update the current singularity file with the new Guppy version  
	a. The current guppy version can be found at [Nanopore Tech Community](https://community.nanoporetech.com/downloads)  
    b. Copy the link next to `Ubuntu 20 GPU`. While SciNet is running under CentOS 7, this container is using Ubuntu 20, as it is what I am most familiar with. There is minimal overhead to running a linux distribution in a container under a linux host ([see this source](https://stackoverflow.com/questions/21889053/what-is-the-runtime-performance-cost-of-a-docker-container))  
	c. Navigate to `pipeline/setup/` and edit the `Singularity` file  
    d. Update the line `DOWNLOAD_LINK=` with the link you have copied  
    e. Save and exit this file by typing `CTRL + x` -> `y` -> `ENTER`  
7. Run `singularity build --remote guppy_new_container.sif pipeline/setup/Singularity`
	a. This will take a bit of time (5 to 10 minutes)
    b. We are not going to overwrite the old file until we are sure the new one is able to build
8. **Assuming no errors occurred**, we will overwrite the old container with the new one
	a. Run `mv guppy_new_container.sif guppy_container.sif`
9. A new guppy container is available.
    a. To test it, see the following section
	

### Running Guppy in a Singularity Container
The guppy container can be ran as an executable, even outside snakemake. Use the following format to interact with the container  
This assumes the guppy container is located at `/project/brookings_minion/guppy_container.sif`
`singularity exec /project/brookings_minion/guppy_container.sif [COMMAND]`  
1. Examples  
   a. Getting help menu of guppy_barcoder: `singularity exec /project/brookings_minion/guppy_container.sif guppy_barcoder --help`  
   b. Running a fast5 through the container: `singularity exec /project/brookings_minion/guppy_container.sif guppy_basecaller -i /project/brookings-minion/my_fast5_file.fast5 -o /project/brookings-minion/my_ouput_folder`
   c. See [singularity exec](https://singularity.lbl.gov/docs-exec), [singularity run](https://singularity.lbl.gov/docs-run), and [singularity shell](https://singularity.lbl.gov/docs-shell) for more information and various methods of interacting with singularity containers
