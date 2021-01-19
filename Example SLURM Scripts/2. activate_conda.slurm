#!/bin/bash
#SBATCH --job-name="Activate Conda"	# The name of your job
#SBATCH -p debug                	# Name of the queue you are submitting to (Community Partitions are: short, medium, long, long60, mem, longmem, mem768, debug)
#SBATCH -N 1                  		# The number of nodes in this job

#SBATCH -n 1						# The number of cores/tasks in this job. Queues max out on number of cores available per node
									# See https://scinet.usda.gov/guide/ceres/#partitions-or-queues for help

#SBATCH -t 00:00:30           		# The amount of time you would like allocated for this job (hours:minutes:seconds)
#SBATCH --mail-user=""   			# Enter your email address if you would like to receive alerts

#SBATCH --mail-type=BEGIN,END,FAIL 	# When would you like to receive an email?
									# Options: BEGIN,END,FAIL,TIME_LIMIT

#SBATCH -o "activate_conda_output"	# standard out %j adds job number to output file name and %N adds the node name
#SBATCH -e "activate_conda_output"	# optional, it prints out standard error
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


# In this script we are simply going to demonstrate the process of activating a conda environment

# First, we must source our .bashrc file.
# This is required because SLURM jobs are sent to a separate server, which does not know about our conda environments.
source ~/.bashrc

# Next, we are going to show all conda environments. Because the SLURM job does not have any active environments at the start, it will show an asterisk (*) on the `base` environment
printf "Before activation\n"
conda info --envs


# Now we are going to activate the mapt_pipeline conda environment
# This may fail for several reasons
# 1) You do not have a conda environment named `mapt_pipeline`
# 2) A conda environment IS named mapt_pipeline, but cannot be found by conda.
# 3) To fix this, replace mapt_pipeline with the `prefix` location of your mapt_pipeline equivalent
conda activate mapt_pipeline


# Show the active environments after activation
printf "\n"
printf "After activation\n"
conda info --envs
