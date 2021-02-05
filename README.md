# MAPT

### Naming
This has been my most notable, and in-depth, project to date. As a result, I feel it should have a proper name.
Due to the nature of this project revolving around mapping noisy microbial reads (hopefully) down to the species level, 
the following name seemed appropriate: Microbial Automated Processing Tool, or MAPT for short.


### Usage

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
    a. `git clone https://github.com/JoshLoecker/pipeline`  
    b. This will create a new folder `pipeline` in your current directory. Move this to the folder you would like it to be contained in  
2. Within the `pipeline` folder, edit the `environment.yaml` file.  
    a. Change the last line to the location where you would like to store your conda environments  
    b. On SciNet, storage space in your home directory is limited. It is [strongly recommended](https://scinet.usda.gov/guide/ceres/#quotas-on-home-and-project-directories) to store environments in your `/project` directory.  
    c. By saving environments in your `/project` directory, it is possible for everyone on your team to use the same environment, even if multiple people are running jobs at the same time.
    d. An example of this is creating a conda environment with the name `pipeline` under the `/project/my-project` directory. Edit the `prefix: ` line as follows  
    1) `prefix: /project/-my-project/pipeline`
   
If this is the first time setting up the conda environment for your group, continue to Step 3  
If step 2 has been completed for your group already, you are ready to activate the conda environment. Proceed to Step 4  

3. Navigate to the `setup` folder where you first downloaded the pipeline (from Step 1)  
    a. Call the `setup.sh` script by running `./setup.sh`  
    b. This will run through setting up the conda environment  
4. To activate the conda environment, the `prefix: ` must be known.
    a. Simply type `conda activate /path/to/prefix/name` to activate the conda environment  
    a. To deactivate the environment, type `conda deactivate`.  
5. To be able to use the Guppy suite of tools, a singularity container is available
   a. This container will be used by snakemake during its run.  
   b. For more information, view Running Guppy in a singularity container at the bottom of this page
6. The final step is to edit several lines within the `pipeline/config.yaml` file  
    a. First, set the `results`, `basecall_files` (or `barcode_files`), `reference_database`, and `guppy_container` to their appropriate locations  
	1) `results` is where you would like results of the pipeline to be stored.  
	2) `basecall_files` holds your fast5/fastq files for Guppy.
	   a) If you are going to skip basecalling and do only barcoding, the `basecall_files` config can be left empty (i.e. `""`)
	3) `reference_database` is the database you will be using with MiniMap for alignments.  
	4) `guppy_container` is the location of the guppy container you will be using.  
	
	b. Set any other values required under the `DEFAULT VALUES` section.

**NOTE**: Any values not changed during in the `pipeline/config.yaml` file will remain as-is during the execution of the pipeline. This **will** lead to unexpected outcomes.

### Running the Pipeline
Once these steps are done, the pipeline is ready to run. The pipeline can be run in several methods

1. Interactive Runs  
    a. If you would like to see the output of jobs as they happen, or you have a short job you would like to ensure is working, Interactive Runs can be useful  
    b. Follow the [guide](https://scinet.usda.gov/guide/ceres/#interactive-mode) here for help on how to set up an interactive run
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
	b. A dry-run allows you to see what steps need to be done, and ensures preliminary configuration is set up correctly.  
    c. To perform a dry-run, activate the `mapt_pipeline` conda environment by running the following command: `conda activate mapt_pipeline`  
    d. Next, call snakemake with a dry-run: `snakemake -j 1 --use-singularity -n`  
    1) `snakemake`: Call snakemake  
	2) `-j 1` (or `--cores 1`): Use 1 core for the dry-run (we do not want to bog-down the login node on the cluster. More cores will not speed this process up)  
	3) `--use-singularity`: Use singularity in the pipeline. This allows us to use the Guppy container  
	4) `-n` (or `--dry-run`): This is the dry-run flag for snakemake  
	
	d. This will output a fair amount of information, showing what rules need to be completed. The pipeline can then be ran following point 1 or 2 above.  

## Notes to Future Maintainers
### How to build a new Guppy container

Building a new guppy container is relatively simple
1. Go to [SyLabs Cloud Builder](https://cloud.sylabs.io/builder)
	a. Sign in with one of the options available
	b. Click on your username in the top right corner
	c. On the drop down, click "Access Tokens"
5. Generate a new Access Token
	a. Copy this token, and save it in a secure location. It will not be able to be accessed again
6. [Log in to SciNet](https://scinet.usda.gov/guide/ceres/#system-access)
7. Navigate to /project/brookings_minion/ on SciNet
	a. There is most likely a file named `guppy_container.sif`
	b. We are going to be updating this file
9. Run `singularity remote login` and follow the instructions to enable the remote build server
10. Update the current singularity file with the new Guppy version
	a. The current guppy version can be found at [Nanopore Tech Community](https://community.nanoporetech.com/downloads)
    b. Update the line `GUPPY_VERSION=` with the new version number
    c. Save and exit this file by typing `CTRL + x` -> `y` -> `ENTER`
10. Run `singularity build --remote guppy_new_container.sif pipeline/setup/Singularity`
	a. This will take a bit of time (5 to 10 minutes)
    b. We are not going to overwrite the old file until we are sure the new one is able to build
11. **Assuming no errors occurred**, we will overwrite the old container with the new one
	a. Run `mv guppy_new_container.sif guppy_container.sif`
12. A new guppy container is available.
    a. To test it, see the following section
	

### Running Guppy in a singularity container
The guppy container can be ran as an executable, even outside snakemake. Use the following format to interact with the container  
`singularity exec [GUPPY_CONTAINER_NAME] [COMMAND]`  
1. Examples  
   a. Getting help menu of guppy_barcoder: `singularity exec guppy_container guppy_barcoder --help`  
   b. Running a fast5 through the container: `singularity exec guppy_container.sif guppy_basecaller -i /project/brookings-minion/my_fast5_file.fast5 -o /project/brookings-minion/my_ouput_folder`
   c. See [singularity exec](https://singularity.lbl.gov/docs-exec), [singularity run](https://singularity.lbl.gov/docs-run), and [singularity shell](https://singularity.lbl.gov/docs-shell) for more information and various methods of interacting with singularity containers
