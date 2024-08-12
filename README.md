# Overview 

1. Installation
2. Files

## Installation
Julia v. >1.9 is required https://julialang.org/
for the Notebooks Jupyter Lab is required https://jupyter.org/

The project uses Packages, that can be used using the Julia package manager. (Opening the Julia REPL `julia` in the directory that contains the package then the Package Manager with `]` and activating the project (with `activate BeyondHulten`)
All the packages will then be installed with the `instantiate` command. The list of packages can be found in `Project.toml`

Also the I-O file has to be inserted into the /data directory

## Files

The `Analysis.ipynb` file contains code from the original B&F paper, that was then translated into Julia. A more condensed replication can be found in `Translation.ipynb`. In the Notebook `DemandShocks.ipynb` the model works on the data from germany and with demand shocks already implemented. This is probably the most interresting to look at. The model can be found in `src` directory.



