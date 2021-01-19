#!/bin/bash
#SBATCH --job-name="Dry-Run"   	# The name of your job
#SBATCH -p debug                	# Name of the queue you are submitting to (Community Partitions are: short, medium, long, long60, mem, longmem, mem768, debug)
#SBATCH -N 1                  		# The number of nodes in this job

#SBATCH -n 1						# The number of cores/tasks in this job. Queues max out on number of cores available per node
									# See https://scinet.usda.gov/guide/ceres/#partitions-or-queues for help

#SBATCH -t 00:00:30           		# The amount of time you would like allocated for this job (hours:minutes:seconds)
#SBATCH --mail-user=""   			# Enter your email address if you would like to receive alerts

#SBATCH --mail-type=BEGIN,END,FAIL 	# When would you like to receive an email?
									# Options: BEGIN,END,FAIL,TIME_LIMIT

#SBATCH -o "snakemake_output"		# standard out %j adds job number to output file name and %N adds the node name
#SBATCH -e "snakemake_output"		# optional, it prints out standard error
									# String replacement options are as follows
										# \\: Do not process any of the replacement symbols.
										# %%: The character "%".
										# %A: master job allocation number.
										# %a: Job array ID (index) number.
										# %J: jobid.stepid of the running job. (e.g. "128.0")
										# %j: jobid of the running job.
										# %N: short hostname. This will create a separate IO file per node.
										# %n: Node identifier relative to current job (e.g. "0" is the first node of the running job) This will create a separate IO file # per node.
										# %s: stepid of the running job.
										# %t: task identifier (rank) relative to current job. This will create a separate IO file per task.
										# %u: User name.
										# %x: Job name.

# NOTE
# This script assumes you have set up your conda environment for a snakemake run. If this is not done properly, snakemake will fail. This script is set up to use my own configuration files.
# If you would like to view these, navigate to `/90daydata/shared/ncarl_minion/joshl/pipeline` and view the `config.yaml` file


# activate conda environment
# This may fail for several reasons
# 1) You do not have a conda environment named `mapt_pipeline`
# 2) A conda environment IS named mapt_pipeline, but cannot be found by conda.
# 3) To fix this, replace mapt_pipeline with the `prefix` location of your mapt_pipeline equivalent
conda activate /project/brookings_minion/conda-envs/mapt_pipeline

# change to pipeline directory
cd /project/brookings_minion/pipeline/

# show snakemake with a dry-run
snakemake -j all --use-singularity -n
printf "\n"
printf "pipeline complete"

# NOTE
# This will not work as-is due to extra setup being required in the following file
# /project/brookings_minion/pipeline/config.yaml
