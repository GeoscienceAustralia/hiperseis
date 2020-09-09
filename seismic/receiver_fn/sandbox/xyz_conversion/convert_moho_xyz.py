"""
Non-workflow script to convert XYZ Moho data from Kennet into 
point datasets compatible with pointsets2grid.py.

XYZ files should be in format:
    STA LAT LON DEPTH WEIGHT
"""
import os

import click
import pandas as pd 

@click.command()
@click.argument('inpath', type=click.Path(exists=True), required=True)
def main(inpath):
    """
    Convert XYZ Moho data to a format ingestible by pointsets2grid.py.

    Output file(s) will be placed in the working directory, with name
    '{filename}_converted.csv'.

    Parameters
    ----------
    inpath: path
        Path to an '.xyz' file to convert, or a directory. If a 
        directory is provided then all '.xyz' files in the directory
        will be converted.
    """
    if os.path.isdir(inpath):
        files = [os.path.join(inpath, f) for f in os.listdir(inpath) 
                 if os.path.isfile(os.path.join(inpath, f)) and f.endswith('.xyz')]
    else:
        files = [inpath]

    for f in files:
        df = pd.read_csv(f, names=['Sta', 'Lat', 'Lon', 'Depth', 'Weight'],
                         delim_whitespace=True, skiprows=1)
        df = df.reindex(columns=['Sta', 'Lon', 'Lat', 'Depth', 'Weight'])
        outname = os.path.splitext(os.path.basename(f))[0] + '_converted.csv'
        with open(outname, 'w') as outfile:
            outfile.write('# Sta,Lon,Lat,Depth,Weight\n')
            df.to_csv(outfile, header=False, index=False)


if __name__ == '__main__':
    main()
