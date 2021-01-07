#!/usr/bin/env bash
# This file will create a conda environment using the environment file in the parent directory
# If the conda environment already exists, nothing will be done

# Get the name of the environment to create
# Used for maintainability, in case the environment name changes nothing will break
# print name of file, redirecting output to grep. search file for 'name: '
env_name=$(tail < "../environment.yaml" | grep "name: ")
env_name="${env_name:6}" # take data after 'name: ' (i.e. environment name)

# Get a list of current environments
current_conda_envs=$(conda env list)

# we have not found an environment with the name 'pipeline', we need to create it
if [[ "$current_conda_envs" != *"$env_name"* ]]; then
    # Create a new conda environment from the "environment.yaml" file in the parent directory
    conda env create --file "../environment.yaml"
# we have found a conda environment with the name 'pipeline', check if it has the correct packages
else
    printf "Found the required conda environment: %s\n" "$env_name"
fi

if [[ "$?" == 0 ]]; then
    printf "\n"
    printf "Setup complete.\n"
    printf "You can list all environments with 'conda info --envs'\n"
    printf "The MAPT pipeline will be using the environment %s, as shown above\n" "$env_name"
fi
