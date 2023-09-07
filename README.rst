HiPerSeis: High Performance Software Package for Seismology Data/Metadata Processing and Analysis
=================================================================================================

|Build Status| |Coverage Status| |Documentation Status|

How to Cite
===========

If you use this software in a scientific publication, we'd very much appreciate if you could cite the following papers:

-  Hassan, R., Hejrani, B., Medlin, A., Gorbatov, A. and Zhang, F., 2020. High-performance seismological tools (HiPerSeis). In: Czarnota, K., Roach, I., Abbott, S., Haynes, M., Kositcin, N., Ray, A. and Slatter, E. (eds.) Exploring for the Future: Extended Abstracts, Geoscience Australia, Canberra, 1–4. https://ecat.ga.gov.au/geonetwork/srv/eng/catalog.search#/metadata/135095
   

Overview
========

- Home Page: https://github.com/GeoscienceAustralia/hiperseis

- Documentation: http://hiperseis.readthedocs.io/en/develop/

- Wiki Pages: https://github.com/GeoscienceAustralia/hiperseis/wiki



Current Contacts
================

- Rakib Hassan: rakib.hassan@ga.gov.au 

- Alexei Gorbatov: alexei.gorbatov@ga.gov.au

- Babak Hejrani: babak.hejrani@ga.gov.au


System Requirements
==========================

- Python 3.6 (recommended)

Setup Guide
=================================

1. First, obtain the source code from `Github repository <https://github.com/GeoscienceAustralia/hiperseis>`_

-  ``git clone https://github.com/GeoscienceAustralia/hiperseis.git``
- ``cd hiperseis``
- ``git submodule update --init --recursive``

2. HiPerSeis does not provide an installation script due to the number of dependencies involved, some of which require low-level libraries to be available on the host machine. Instead, shell scripts are provided in ``hiperseis/setup_scripts`` for Linux, OSX and Windows for installing dependencies through a combination of Conda and Pip. A shell script is provided for NCI GADI, tailored exclusively for the current list of low-level HPC libraries e.g. MPI, HDF5, etc. available on the system.

3. To use HiPerSeis in the checked out location, you will need to add the root HiPerSeis folder to your PYTHONPATH variable. For example, if you checked out HiPerSeis to the folder `dev/hiperseis` relative to your home directory, then in a `bash` shell you need to execute the following shell command: ``export PYTHONPATH=$HOME/dev/hiperseis``.  This needs to be done for each command shell session, or added to ``.bashrc`` or its equivalent.

Third Party Library Dependencies
================================

Certain modules require specific third party (non-Python) libraries to be installed
on the host system. For example, scripts that convert to sc3ml format also require Seiscomp3 to be
installed and to be visible in the PATH.

License
===============

HiPerSeis is licensed under the GPL version 3



.. |Build Status| image:: https://travis-ci.org/GeoscienceAustralia/hiperseis.svg?branch=develop
   :target: https://travis-ci.org/GeoscienceAustralia/hiperseis
   
.. |Coverage Status| image:: https://coveralls.io/repos/github/GeoscienceAustralia/hiperseis/badge.svg
   :target: https://coveralls.io/github/GeoscienceAustralia/hiperseis

.. |Documentation Status| image:: https://readthedocs.org/projects/hiperseis/badge/?version=develop
   :target: http://hiperseis.readthedocs.io/en/develop/

