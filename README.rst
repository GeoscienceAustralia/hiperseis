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

- Linux OS
- Python 2.7
- Python 3.5 or higher


Setup Guide for Developers
==========================

1. Install Python environment and dependency packages. See Wiki Pages: https://github.com/GeoscienceAustralia/hiperseis/wiki

2. Obtain the source code from https://github.com/GeoscienceAustralia/hiperseis

-  ``git clone https://github.com/GeoscienceAustralia/hiperseis.git``
- ``cd hiperseis``
- ``git submodule init``
- ``git submodule update``

   - ``pip install -v --user -e .`` (into user's own home ~/.local/lib/python2.7/site-packages/)


If you are using the library without ``pip install``-ing it, then make sure the ``hiperseis`` folder
is in your ``PYTHONPATH`` environment variable.  E.g. in bash shell, if ``~/dev/hiperseis`` is where
``hiperseis`` was checked out:

- ``export PYTHONPATH=~/dev/hiperseis:$PYTHONPATH``


License
===============

HiPerSeis is licensed under the GPL version 3



.. |Build Status| image:: https://travis-ci.org/GeoscienceAustralia/hiperseis.svg?branch=develop
   :target: https://travis-ci.org/GeoscienceAustralia/hiperseis
   
.. |Coverage Status| image:: https://coveralls.io/repos/github/GeoscienceAustralia/hiperseis/badge.svg
   :target: https://coveralls.io/github/GeoscienceAustralia/hiperseis

.. |Documentation Status| image:: https://readthedocs.org/projects/hiperseis/badge/?version=develop
   :target: http://hiperseis.readthedocs.io/en/develop/

