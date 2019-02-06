# ateam-tools
Shared modeling and analysis tools for internal use by the Anastassiou research team

For active development, install via `pip install -e <path>` or `setup.py develop`.

For much of the functionality to work, you will need the ateam branch of BMTK, 
at https://github.com/tmchartrand/bmtk.git.

# Package organization and directory structure
```
ateam-tools
├── ateam: base directory for python code files, organized by function
|    ├── analysis
|    ├── sim
|    ├── ...etc 
├── examples: scripts or notebooks that demonstrate use of our tools, organized by application
|   ├── subfolders organized by project or topic
├── scripts: tools for running analyses from the command line
```
Some organizing principles:
- Don't duplicate or copy and paste code! (Make a function instead.)
- Code that has distinct functions (sim vs analysis) or distinct sets of dependencies (BMTK vs aiephys analysis) should be grouped in separate folders. (These distinctions tend to go together.)
- Try to minimize new folders (packages), but generally keep your code in your own modules (.py files)
- Try to minimize new dependencies if possible. If you must use a new package, try to separate code that needs it from code that doesn't into different modules.
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

The environment also contains a full jupyter notebook server, so the `jupyter notebook` command
will start a server running the shared environment.
Alternatively, you may want to use a jupyter plugin to allow jupyter to access and switch between multiple conda environments:

nb_conda_kernels (https://docs.anaconda.com/anaconda/user-guide/tasks/use-jupyter-notebook-extensions/)

or 
jupyter_environment_kernels 
(https://github.com/Cadair/jupyter_environment_kernels)