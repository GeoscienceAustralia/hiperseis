#!/usr/bin/env python

import os
import glob
import subprocess

contact_email = 'andrew.medlin@ga.gov.au'

for f in glob.glob('*.dat'):
    fname, ext = os.path.splitext(f)
    dname = fname + "_OUT"
    if not os.path.isdir(dname):
        print('os.mkdir({})'.format(dname))
        os.mkdir(dname)
        cmd = ['qsub', '-M', contact_email, '-N', fname,
               '-v', 'INFILE={},OUT={}'.format(f, dname), './run_rf.sh']
        print(' '.join(cmd))
#        subprocess.run(['echo'] + cmd, check=True)
        subprocess.run(cmd, check=True)
    else:
        print('Skipping {}'.format(fname))
    # end if
# end for

