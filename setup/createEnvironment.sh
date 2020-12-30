#!/usr/bin/env bash
# This file will create a conda environment using the environment file in the parent directory
# If the conda environment already exists, nothing will be done

project_base_directory="/90daydata/shared/ncarl_minion/joshl/"


# try to load miniconda
load_miniconda_error_status=$(module load miniconda 2>&1)
# we are on SciNet, able to load module miniconda
if [[ $load_miniconda_error_status == 0 ]]; then
    printf "Loaded miniconda module"

    # try to initialize conda
    conda_init_status=$(conda init 2>&1)
    if [[ "$conda_init_status" != 0 ]]; then
        printf "You need to run the command 'conda init' before continuing\n"
        exit 1
    elif [ "$conda_init_status" == 0 ]; then
        printf "Conda is initialized"
    fi
fi

# Get the name of the environment to create
# Used for maintainability, in case the environment name changes nothing will break
# print name of file, redirecting output to grep. search file for 'name: '
env_name=$(< "../environment.yaml" tr ' ' _ | grep "name:")
env_name="${env_name:6}"  # take data after 'name: ' (i.e. environment name)

# Get a list of current environments
current_conda_envs=$(conda env list)

# we have not found an environment with the name 'pipeline', we need to create it
if [[ "$current_conda_envs" != *"$env_name"* ]]; then
    # Create a new conda environment from the "environment.yaml" file in the parent directory
    conda env create --file "../environment.yaml" --prefix "$project_base_directory"
# we have found a conda environment with the name 'pipeline', check if it has the correct packages
else
    printf "Found the required conda environment: %s\n" "$env_name"
fi

printf "\n"
printf "Setup complete.\n"
printf "You can list all environments with 'conda info --envs'\n"
printf "You may need to run 'module load miniconda' first, if running on a cluster\n"
printf "The MAPT pipeline will be using the environment %s, as shown above\n" "$env_name"

