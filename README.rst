===============
Passive-Seismic
===============

Software and metadata for GA passive seismic project
----------------------------------------------------

.. image:: https://circleci.com/gh/GeoscienceAustralia/passive-seismic.svg?style=shield
    :target: https://circleci.com/gh/GeoscienceAustralia/passive-seismic


Dependencies
------------
TODO


Install ``passive-seismic``
---------------------------
TODO


Tests
-----

Tests require ``ELLIPCORR`` env variable to point to the ``ellip-corr``
directory:

.. code:: console

    $ cd /path/to/passive-seismic
    $ export ELLIPCORR=$PWD/ellip-corr

Then run tests without `coverage`:

.. code:: console

    $ pytest tests
    $ # or using the makefile
    $ make test

Or with `coverage`:

.. code:: console

    $ make coverage
