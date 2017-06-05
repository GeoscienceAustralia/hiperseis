Creating a Python virtualenv on ANTELOPE for exporting events XML
=================================================================

This is a quick guide to exporting the GA ANTELOPE database events
information into QuakeML. These instructions are tailored to the GA ANTILOPE
dev system and may need further tuning on the prod system.

These instructions assume you are using bash shell.

----------------
Pre-installation
----------------

1. Install ``virtualenv``:

   .. code:: bash

       $ pip install  --user virtualenv

   This will install ``virtualenv`` binary in ``~/.local/bin/virtualenv``.


2. Clone the ``passive-seismic`` repository into your home directory, or
another directory of your choice ( I chose /export/

   .. code:: bash

       $ cd ~
       $ git clone https://github.com/GeoscienceAustralia/passive-seismic.git

2. Unload the icc compiler and default openmpi from the terminal:

   .. code:: bash

       $ module unload intel-cc
       $ module unload intel-fc
       $ module unload openmpi

3. Load the modules required for installation and running:

   .. code:: bash

       $ module load python3/3.4.3 python3/3.4.3-matplotlib
       $ module load hdf5/1.8.14p openmpi/1.8 mpi4py/2.0.0

   (Alternatively, you may wish to add the above lines to your
   ``~/.profile`` file)

4. Now add the following lines to the end of your ``~/.profile`` file:

   .. code:: bash

       export PATH=$HOME/.local/bin:$PATH
       export PYTHONPATH=$HOME/.local/lib/python3.4/site-packages:$PYTHONPATH
       export VIRTUALENVWRAPPER_PYTHON=/apps/python3/3.4.3/bin/python3
       export LC_ALL=en_AU.UTF-8
       export LANG=en_AU.UTF-8
       source $HOME/.local/bin/virtualenvwrapper.sh

5. Install virtualenv and ``virtualenvwrapper`` in a terminal:

   .. code:: bash

       $ pip3 install  --user virtualenv virtualenvwrapper

6. Refresh your environment by sourcing your ``~/.profile`` file:

   .. code:: bash

       $ source ~/.profile

------------
Installation
------------

1. Create a new virtualenv for ``passive-seismic``:

   .. code:: bash

       $ mkvirtualenv --system-site-packages seismic

2. Make sure the virtualenv is activated:

   .. code:: bash

       $ workon seismic

3. Clone ``h5py`` from ``https://github.com/basaks/h5py.git``:

   .. code:: bash

       $ cd ~
       $ git clone https://github.com/basaks/h5py.git
       $ cd ~/h5py
       $ export CC=mpicc
       $ python setup.py configure --mpi --hdf5=/apps/hdf5/1.8.14p
       $ python setup.py install


4. Install ``passive-seismic``:

   .. code:: bash

       $ cd ~/passive-seismic
       $ python setup.py install

5. Once installation has completed, you can run the tests to verify
   everything has gone correctly:

   .. code:: bash

       $ pip install pytest
       $ pytest tests/
