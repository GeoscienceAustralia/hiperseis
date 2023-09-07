# Anaconda-based Setups

1. Setup scripts provided for Linux, OSX and Windows (conda_env_*.sh) are based on
[Anaconda3-2021.11](https://repo.anaconda.com/archive/).
Refer to Anaconda installation instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)

2. After installing Anaconda (v 2021.11), launch the setup shell script for your OS, e.g. for Linux:

 ``source hiperseis/setup_scripts/conda_env_linux.sh env_name || true``

 where, ``env_name`` is the name of the environment to be created (e.g. hiperseis_env) and
 `` || true`` ensures the shell terminal the command is run from does not terminate
 if errors are encountered.

3. Activate the newly created environment:

``conda activate env_name``

Note that the setup script for Windows is a whittled down version compared to the Linux and OSX versions. While most of the dependecies 
are covered, some that require code compilation are excluded, specifically to cater to institutional Windows machines with strict security controlls in place. ``conda_env_windows.sh`` should be launched from Git BASH -- refer to its installation instructions [here](https://gitforwindows.org/).

# GADI

1. The setup script for GADI is launched as:

``sh hiperseis/setup_scripts/setup_env_gadi.sh env_name``

where, ``env_name`` is the name of the environment to be created (e.g. hiperseis_env).

2. Activate the newly created environment:

``source env_name/bin/activate``
