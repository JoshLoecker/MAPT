#!/bin/bash
#SBATCH --job-name="Simple Job"   	# The name of your job
#SBATCH -p debug                	# Name of the queue you are submitting to (Community Partitions are: short, medium, long, long60, mem, longmem, mem768, debug)
#SBATCH -N 1                  		# The number of nodes in this job

#SBATCH -n 1						# The number of cores/tasks in this job. Queues max out on number of cores available per node
									# See https://scinet.usda.gov/guide/ceres/#partitions-or-queues for help

#SBATCH -t 00:00:30           		# The amount of time you would like allocated for this job (hours:minutes:seconds)
#SBATCH --mail-user=""   			# Enter your email address if you would like to receive alerts

#SBATCH --mail-type=BEGIN,END,FAIL 	# When would you like to receive an email?
									# Options: BEGIN,END,FAIL,TIME_LIMIT

#SBATCH -o "simple_job_output"		# standard out %j adds job number to output file name and %N adds the node name
#SBATCH -e "simple_job_output"		# optional, it prints out standard error
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


# This is a very simple SLURM job. We are going to show the current date/time.

# Looking at the #SBATCH above, we are requesting 1 Debug node, with 1 core. This job will run for a maximum of 30 seconds, or until the following commands are finished, whichever comes first.
# We will not receive an email when this is ran, but this can be configured if you would like
# If an email is supplied, we will receive an email at the beginning, end, and any possible failures of the job

# Any output is going to be saved into a file with the prefix `output.####`, where `####` is the job ID

# Start this script by calling `sbatch simple_job.slurm`


date
