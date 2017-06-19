Using the ANU Logfile Script
============================

The ``passive-seismic/convert_logs/decode_datfile.py`` can be used to convert
the binary `.dat` ANU log files into ``jsons``.

Can only use python2.7 for now.

------------------
Python Environment
------------------

Environment has these python packages:

   .. code:: bash

    click (6.7)
    numpy (1.13.0)
    pip (9.0.1)
    setuptools (36.0.1)
    wheel (0.29.0)
    obspy (1.0.3)

------------------
Running the script
------------------

There are two ways to run the ``anulog`` conversions:

#. directly run the python script
#. after ``passive-seismic`` is installed, you can use the command ``anulog``

Checkout the help string:

   .. code:: bash

    $ python decode_datfile.py --help

After ``passive-seismic`` is installed, you can simply use the ``anulog`` command

   .. code:: bash

    $ anulog --help

Both will produce the following help string:

::

 Usage: anulog [OPTIONS] DATFILE

  Program to display contents of the logfile <datfile>.dat

 Options:
  -b, --bad_gps BOOLEAN       Print bad gps info
  -u, --gps_update BOOLEAN    Print the gps update info
  -i, --id_str BOOLEAN        Print bad id strings
  -t, --temperature BOOLEAN   Print bad temperature info
  -a, --all_print BOOLEAN     Print all
  -y, --year INTEGER RANGE    Gpsyear. max(Gpsyear - year) == 1
  -o, --output_dir DIRECTORY  Output dir name. If no output dir is provided,
                              input dir will be used.
  -v, --verbosity [DEBUG|INFO|WARNING|ERROR]
                                  Level of logging
  --help                      Show this message and exit.

A typical log file conversion command is just the following:

   .. code:: bash

    $ python decode_datfile.py datfile
    $ # or
    $ anulog datfile

This will output ``json``s  corresponding to the ``.dat`` ``anulog`` files.
The input `datfile` could be one `.dat` log file, or can be a directory of
``.dat`` logfiles.

---------------------------
Using inside custom scripts
---------------------------

Once ``passive-seismic`` is installed, you have access to the decode
functionality ``decode_anulog`` to use in your script. Import it in your script
like the following:

   .. code:: bash

    $ In [1]: from convert_logs.decode_datfile import decode_anulog


The function ``decode_anulog`` will output a python ``dict`` corresponding to
 the binary ``anulog`` file.

-------------------
Parallel conversion
-------------------

The ``decode_datfile.py``/``anulog`` code uses ``multiprocessing``. To use
multiprocessing in your script you can use the following:

   .. code:: bash

    $ In [2]: from joblib import Parallel, delayed
    $ In [3]: datfiles = glob.glob(os.path.join(datfile_dir, '*.dat'))
    $ In [4]: log_dicts = Parallel(n_jobs=-1)(delayed(decode_anulog)(
                  d, bad_gps, id_str, gps_update, temperature, all_print, year)
                                          for d in datfiles)
