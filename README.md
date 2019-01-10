# ateam-tools
Shared modeling and analysis tools for internal use by the Anastassiou research team

For active development, install via `pip install -e <path>` or `setup.py develop`.

For much of the functionality to work, you will need the ateam branch of BMTK, 
at https://github.com/tmchartrand/bmtk.git.

# Shared Conda environment for BMTK and NEURON
In order to have a stable, reproducible platform for testing and running simulations,
we have created a standalone installation of the ateam-tools, BMTK, and NEURON software stack.
The installation is packaged as a conda environment, and automatically configures the necessary 
paths for running NEURON when activated.

## Setup of bmtk_ateam Conda environment
Simply run the shell script `/allen/aibs/mat/ateam_shared/software/add_shared_conda.sh` to add this 
shared conda environment to your local install of conda/miniconda.

## Using the bmtk_ateam Conda environment
Following setup, simply call `conda activate bmtk_ateam` or `source activate bmtk_ateam` to 
activate the environment and use.

Although this environment does not contain a full jupyter notebook server, it can also be 
used for jupyter notebooks. To do so requires a plugin to allow jupyter to access conda environments:

nb_conda_kernels (https://docs.anaconda.com/anaconda/user-guide/tasks/use-jupyter-notebook-extensions/)

or 
jupyter_environment_kernels 
(https://github.com/Cadair/jupyter_environment_kernels)