#!/bin/bash
#SBATCH --job-name="TMPDIR"			# The name of your job
#SBATCH -p debug                	# Name of the queue you are submitting to (Community Partitions are: short, medium, long, long60, mem, longmem, mem768, debug)
#SBATCH -N 1                  		# The number of nodes in this job

#SBATCH -n 1						# The number of cores/tasks in this job. Queues max out on number of cores available per node
									# See https://scinet.usda.gov/guide/ceres/#partitions-or-queues for help

#SBATCH -t 00:00:30           		# The amount of time you would like allocated for this job (hours:minutes:seconds)
#SBATCH --mail-user=""   			# Enter your email address if you would like to receive alerts

#SBATCH --mail-type=BEGIN,END,FAIL 	# When would you like to receive an email?
									# Options: BEGIN,END,FAIL,TIME_LIMIT

#SBATCH -o "tmpdir_README"        		# standard out %j adds job number to output file name and %N adds the node name
#SBATCH -e "tmpdir_README"        		# optional, it prints out standard error
									# String replacement options are as follows
										# \\: Do not process any of the replacement symbols.
										# %%: The character "%".
										# %A: master job allocation number.
										# %a: Job array ID (index) number.
										# %J: jobid.stepid of the running job. (e.g. "128.0")
										# %j: jobid of the running job.
										# %N: short hostname. This will create a separate IO file per node.
										# %n: Node identifier relative to current job (e.g. "0" is the first node of the running job) This will create a separate IO 		file # per node.
										# %s: stepid of the running job.
										# %t: task identifier (rank) relative to current job. This will create a separate IO file per task.
										# %u: User name.
										# %x: Job name.

cat << _EOF
According to the following SciNet website, it is wise to copy large amounts of data to the ~9TB SSD storage for large jobs
These SSDs are faster than the hard drives used in the ~/project/ folder, and will result in a faster computation time overall
https://scinet.usda.gov/guide/ceres/#local-scratch-space-on-large-memory-nodes

This guide will attempt to show the differences between the current file system and that available under TMPDIR
Testing shows the following results:
/90daydata/:    175 MB/s write speed (9 minutes 31 seconds for 100 GB of data
\$TMPDIR:       195 MB/s write speed (8 minutes 32 seconds for 100 GB of data)

While at first glance this is not a significant increase, the benefit also comes from higher input/output operations being available on an SSD.
_EOF

printf "\n\n"
printf "Print the path to our temporary directory (known as TMPDIR in the future)\n"
echo "$TMPDIR"

printf "\n\n"
printf "View the contents of TMPDIR\n"
ls "$TMPDIR"

printf "\n\n"
printf "Copying files to the temporary directory just like any other operation\n"
printf "cp 4_copy_to_tmpdir.slurm %s\n" "$TMPDIR"
cp 4_copy_to_tmpdir.slurm $TMPDIR

printf "\n\n"
printf "View files in the temp directory after we have copied a file\n"
ls "$TMPDIR"

printf "\n\n"
cat << _EOF
If we make a change to a file in the $TMPDIR, we must copy all files back from it
The files we copy are not located on the current filesystem - once this script exits, the data will be lost
We are going to create a file in TMPDIR to show the separation of file systems
_EOF
tmpdir_file="this_file_is_in_tmpdir"
touch "$tmpdir_file"

printf "\n\n"
printf "This is currently inside TMPDIR\n"
ls "$TMPDIR"

printf "\n\nAs we can see, a file named %s is within" "$tmpdir_file\n"
printf "However, this file is not in the current /project/ directory. Watch below\n"
printf "Contents of %s" "$(pwd)"
ls "$(pwd)"

printf "\n\n"
printf "This isn't a problem. As stated previously, we are able to copy files between the TMPDIR and our current directory\n"
printf "We will now copy the %s file to our current directory (%s)" "$tmpdir_file" "$(pwd)"
cp "$TMPDIR/$tmpdir_file" "$(pwd)"

printf "\n\n"
printf "We can now see there is a file named %s in our current directory\n" "$tmpdir_file"
printf "ls \$pwd\n"
ls "$(pwd)"


cat << _EOF
However, just because a file was copied once does not mean changes to the file will be synced.
If we type a string into the file '$tmpdir_file', those changes will not propagate into the current file
We are going to do the following: 'cat "This string is in \$tmpdir_file" > \$tmpdir_file'
The string 'This string is in $tmpdir_file' will be found within the file
_EOF
cat "This string is in the \$tmpdir_file" > "$tmpdir_file"

printf "\n\n"
printf "We can see the changes by viewing the file, as such:\n"
cat $tmpdir_file

printf "\n\n"
printf "Viewing the contents of the file in our current directory, we can see those changes are not present:\n"
cat "'$(pwd)/$tmpdir_file'"

printf "\n\n"
printf "This shows that files (and changes to the data) must be copied TO the tmpdir at the beginning of slurm scripts and FROM the tmpdir end of slurm scripts\n"
printf "The temporary file that was placed in this directory to show the differences in file systems has been removed\n\n"
printf "If this file was viewed with 'cat', please scroll to the top\n. Thank you"

rm "$(pwd)/$tmpdir_file"


