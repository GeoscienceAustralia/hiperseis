HiPerSeis: High Performance Software Package for Seismology Data/Metadata Processing and Analysis
=================================================================================================

|Build Status| |Coverage Status| |Documentation Status|


Overview
========

- Home Page: https://github.com/GeoscienceAustralia/hiperseis

- Documentation: http://hiperseis.readthedocs.io/en/develop/

- Wiki Pages: https://github.com/GeoscienceAustralia/hiperseis/wiki



Contacts
==========

- Fei Zhang: fei.zhang@ga.gov.au

- Rakib Hassan: rakib.hassan@ga.gov.au

- Andrew Medlin: andrew.medlin@ga.gov.au

- Alexei Gorbatov: alexei.gorbatov@ga.gov.au

- Babak Hejrani: babak.hejrani@ga.gov.au


System Requirements
==========================

- Linux OS, including Ubuntu and CentOS
- Python 3.6 or higher (recommended)
- Python 2.7 (deprected, no longer supported)

Third Party Library Dependencies
================================

Certain modules require specific third party (non-Python) libraries to be installed
on the host system. For example, scripts that convert to sc3ml format also require Seiscomp3 to be
installed and to be visible in the PATH. In most cases, Python libraries that depend on third party
libraries will indicate their dependencies either during attempted installation, or when used at
runtime. Note that the following list includes indirect dependencies, i.e. dependencies that come
from Python libraries used by HiPerSeis.

Current third party dependencies (actual requirements may vary by platform or Python distribution):

- `HDF5 <http://hdfgroup.org/>`_
- MPI, for example `Open MPI <https://www.open-mpi.org/>`_
- `PROJ <https://proj.org/>`_
- `GEOS <https://trac.osgeo.org/geos>`_


Installation Guide for Developers
=================================

1. First, obtain the source code from `Github repository <https://github.com/GeoscienceAustralia/hiperseis>`_

-  ``git clone https://github.com/GeoscienceAustralia/hiperseis.git``
- ``cd hiperseis``
- ``git submodule init``
- ``git submodule update``

2. Install Python environment and dependency packages. See `Wiki Pages <https://github.com/GeoscienceAustralia/hiperseis/wiki>`_

3. To use HiPerSeis in the checked out location, you will need to add the root HiPerSeis folder to your PYTHONPATH variable. For example, if you checked out HiPerSeis to folder `dev/hiperseis` relative to your home directory, then in `bash` shell you need to execute the following shell command: ``export PYTHONPATH=$HOME/dev/hiperseis``.  This needs to be done for each command shell session.

License
===============

HiPerSeis is licensed under the GPL version 3



.. |Build Status| image:: https://travis-ci.org/GeoscienceAustralia/hiperseis.svg?branch=develop
   :target: https://travis-ci.org/GeoscienceAustralia/hiperseis
   
.. |Coverage Status| image:: https://coveralls.io/repos/github/GeoscienceAustralia/hiperseis/badge.svg
   :target: https://coveralls.io/github/GeoscienceAustralia/hiperseis

.. |Documentation Status| image:: https://readthedocs.org/projects/hiperseis/badge/?version=develop
   :target: http://hiperseis.readthedocs.io/en/develop/

