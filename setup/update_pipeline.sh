# This is a simple script used to update the pipeline
# It will do the following
#   1) Move the user's configuration file to the home directory
#       a) This is done so git does not complain about committing file changes
#   2) Pull the latest changes from GitHub
#   3) Move the configuration file back into the pipeline directory

cd ../
snakefile_present=$(ls snakefile > /dev/null 2>&1)

if [[ "$snakefile_present" != "0" ]]; then
    echo "Execute this script inside the 'pipeline/setup' directory that was cloned from GitHub"
    exit 1
fi

mv config.yaml ~/.config.yaml
git pull
mv ~/.config.yaml ./config.yaml
