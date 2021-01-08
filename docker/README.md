The Dockerfile was moved to its own folder so it did not attempt to pull in all folders to it's build context.
This drastically reduces the number of items listed in the .dockerignore (i.e. 0 items are now listed)
