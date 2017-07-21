Creating a Python virtualenv on ANTELOPE for exporting events XML
=================================================================

This is a quick guide to exporting the GA ANTELOPE database events
information into QuakeML. These instructions are tailored to the GA ANTILOPE
dev system and may need further tuning on the prod system.

These instructions assume you are using bash shell.

-----------------
Create Virtualenv
-----------------

1. Install ``virtualenv``:

   .. code:: bash

       $ pip install  --user virtualenv

   This will install ``virtualenv`` binary in ``~/.local/bin/virtualenv``.


2. Create and activate your ``virtualenv``. I chose
to install my ``virtualenv``s in ``/export/development/sudipta/venvs``
due to low space allocation in my home directory.

   .. code:: bash

       $ cd /export/development/sudipta/
       $ mkdir venvs
       $ ~/.local/bin/virtualenv --system-site-packages /export/development/sudipta/venvs/antelope
       $ source /export/development/sudipta/venvs/antelope/bin/activate

3. Upgrade ``obspy`` as the ``obspy`` in the ANTELOPE system could be very old.
Also install ``lxml`` without using the binaries.

   .. code:: bash

       $ pip install git+https://github.com/basaks/obspy.git -U --no-deps
       $ pip install lxml --no-binary :all:


-----------------
Events to QuakeML
-----------------

1. Clone the ``antelope_contrib`` repository into your home directory, or
another directory of your choice (I chose ``/export/development/sudipta/``
due to low space allocation in my home directory.)

   .. code:: bash

       $ cd /export/development/sudipta/
       $ git clone https://github.com/GeoscienceAustralia/antelope_contrib.git

2. Test one event can be converted to QuakeML

   .. code:: bash

       $ cd antelope_contrib/bin/export/events
       $ python ga_event2qml.py -s schemas/QuakeML-BED-1.2.rng -o event_id.xml db_loc event_id -d

Note the ``db_loc`` has to point to the correct db, and ``event_id`` has
to exist in the antelope database. The ``ga_event2qml.py`` converts the
event corresponding to the ``event_id`` into the ``QuakeML`` ``event_id
.xml``. The script also runs a validation routine on the generate
``QuakeML``.


----------------------
QuakeML to Seiscomp3ML
----------------------
Once the single event test is successful, proceed to covert all of the events
 in the ANTELOPE database into Seiscomp3 compatible XML.

Convert all events in the ANTELOPE database:

   .. code:: bash

       $ python extract_events.py -s schemas/QuakeML-BED-1.2.rng -o sc3.xml db_path


This will generate the QuakeML files inside the ``outdir`` and a
corresponding ``seiscomp3`` xml file ``sc3.xml``.


------------------
Ingest Seiscomp3ML
------------------

This ``sc3.xml`` can be imported into ``seiscomp3`` using the following command

   .. code:: bash

      $ scdb -i sc3.xml -d mysql://sysop:sysop@localhost/seiscomp3
