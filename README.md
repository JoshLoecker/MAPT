Pipeline (MAPT)
---------------

Naming
------
I would like to name this project as it has been my first, and most in-depth, project to date. 
Due to the nature of this project revolving around mapping noisy microbial reads (hopefully) down to the species level, 
the following name seemed appropriate: Microbial Automated Processing Tool, or MAPT for short.


Usage
-----

This is a workflow/pipeline aimed to streamline post-processing of noisy reads from [Oxford Nanopore MinION Sequencing Data](https://www.nanoporetech.com). The following software will need to be installed in a conda environment. A startup-script is provided in the GitHub repository. This is explained further on this page.

The following software will be used (not necessarily in this order):
- Guppy Basecalling
- Guppy Barcoding
- Cutadapt
- NanoFilt
- isONclust
- spoa
- NanoPlot
- Plotly Graphs

The page to the GitHub and Docker repositories are as follows, if they are needed.
- GitHub: [https://github.com/JoshLoecker/pipeline](https://github.com/JoshLoecker/pipeline)
- Docker: [https://hub.docker.com/repository/docker/joshloecker/pipeline](https://hub.docker.com/repository/docker/joshloecker/pipeline)

Installation & Setup
--------------------

This project was built and testing on the following versions of Singularity and Docker. Reproducibility and stability cannot be guaranteed on earlier versions.
1. Singularity (version 3.6.4-1.el7 or higher)
2. Docker (version 20.10.0 or higher)

[Download singularity](https://singularity.lbl.gov/)  
[Download docker](https://www.docker.com/products/docker-desktop)

Git and Miniconda are required to download the pipeline and create a new Conda environment with the required software

[Download git](https://git-scm.com/downloads)
[Download miniconda](https://docs.conda.io/en/latest/miniconda.html)

1. Start by cloning the GitHub repo  
    a. `git clone https://github.com/JoshLoecker/pipeline`
    b. This will create a new folder `pipeline` in your current directory. Move this to the folder you would like it to be contained in
2. Within the `pipeline` folder, edit the `environment.yaml` file.  
    a. Change the last line to the location where you would like to store your conda environments  
    b. On SciNet, storage space in your home directory is limited. It is [strongly recommended](https://scinet.usda.gov/guide/ceres/#quotas-on-home-and-project-directories) to store environments in your `/project` directory.  
    c. By saving environments in your `/project` directory, it is possible for everyone on your team to use the same environment, even if multiple people are running jobs at the same time.
    d. You may also change the name of the pipeline by changing `name: mapt_pipeline` to `name: MY_NAME` if you choose
   
If this is the first time setting up the conda environment for your group, continue to Step 3
If step 2 has been completed for your group already, you are ready to activate the conda environment. Proceed to Step 4

3. Navigate to the `setup` folder where you first downloaded the pipeline (from Step 1)  
    a. Call the `setup.sh` script by running `./setup.sh`  
    b. This will guide you through a simple setup to install the conda environment  
4. If you know where the `prefix: ` to the conda environment was set, simply type `source activate /path/to/conda/env/name`  
    a. To deactivate the environment, type `conda deactivate`.  
5. To be able to use the Guppy suite of tools, the 
5. The final step is to edit several lines within the `pipeline/config.yaml` file  
    a. First, set the `results`, `basecall_files` (or `barcode_files`), `reference_database`, and `guppy_container` to their appropriate locations  
	1) `results` is where you would like results of the pipeline to be stored.  
	2) `basecall_files` holds your fast5/fastq files for Guppy.
	   a) If you are going to skip basecalling and do only barcoding, the `basecall_files` config can be left empty (i.e. `""`)
	3) `reference_database` is the database you will be using with MiniMap for alignments.  
	4) `guppy_container` is the location of the guppy container you will be using.  
	
	b. Then set any other values required under the `DEFAULT VALUES` section. If these are not changed, they will remain as-is during the pipeline run  

Running the Pipeline
--------------------
Once these steps are done, the pipeline is ready to run. The pipeline can be run in several methods

1. Interactive Runs  
    a. If you would like to see the output of jobs as they happen, or you have a short job you would like to ensure is working, Interactive Runs can be useful  
    b. Follow the [guide](https://scinet.usda.gov/guide/ceres/#interactive-mode) here for help on how to set up an interactive run  
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
	a. A dry-run allows you to see what steps need to be done, and ensures preliminary configuration is set up correctly.  
    b. To perform a dry-run, activate the `mapt_pipeline` conda environment by running the following command: `conda activate mapt_pipeline`  
    c. Next, call snakemake with a dry-run: `snakemake -j 1 --use-singularity -n`  
    1) `snakemake`: Call snakemake  
	2) `-j all` (or `--cores 1`): Use 1 core for the dry-run  
	3) `--use-singularity`: Use singularity in the pipeline. This allows us to use the Guppy container  
	4) `-n` (or `--dry-run`): This is the dry-run flag for snakemake  
	
	d. This will output a fair amount of information, showing what rules need to be completed  

Notes to Future Maintainers
---------------------------
### How to build a new Guppy container
This can be done in one of two ways.
1. Automated Method  
    a. This is probably the most preferred method.  
    b. To start, update the `ARG version=. . .` in the `Dockerfile` to the newest version of guppy. The newest version can be found on [Oxford Nanopore Downloads](https://community.nanoporetech.com/downloads)  
    c. Next, execute the script `buildGuppy.sh` under the `setup` folder  
    d. This will go through the process of downloading the specified version of Guppy in a new docker image, and create a new singularity container from the docker image
    e. Once the script is done, the resulting singularity image needs to be uploaded to SciNet for use on the cluster
2. Manual Method  
	a. This method is more intensive, but may be required if something breaks in the automated method  
    b. First, singularity and docker must be installed on the same machine to update guppy  
    c. Building a new Guppy singularity container is first done by creating a docker image  
    1) `docker build --tag [YOUR TAG]`  
    2) This was done simply because I was most familiar with docker containers    
    3) It may be smart to move the Dockerfile in the `pipeline` repository to a Singularity file  
    4) The `Dockerfile` located in this repository is what the image should be built upon  
	
	d. The singularity container is built in the following manner  
	1) `singularity build --sandbox [GUPPY_CONTAINER_NAME] docker-daemon://[YOUR TAG FROM STEP 2]`  
	2) The `docker-daemon` is used for a local docker image. Local images are generally preferred. This means we do not have to upload the resulting container to Dockerhub, then download it to our local machine  
	3) `singularity build --sandbox [GUPPY_CONTAINER_NAME] docker://[YOUR TAG]` will download a docker container from Dockerhub, if this is preferred.  
	
	e. This will ultimately generate a singularity container with the name `[GUPPY_CONTAINER_NAME]`.  
	f. It can be run as an executable in the following format `singularity exec [GUPPY_CONTAINER_NAME] [COMMAND]`
    1) See [singularity exec](https://singularity.lbl.gov/docs-exec), [singularity run](https://singularity.lbl.gov/docs-run), and [singularity shell](https://singularity.lbl.gov/docs-shell) for more information
