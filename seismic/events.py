# import seismic.cluster.cluster as cluster
import os
import sys
import fnmatch


def recursive_glob(adir, ext='*.xml'):
    """
    source: https://stackoverflow.com/a/2186565/3321542
    """
    matches = []
    for root, dirnames, filenames in os.walk(adir):
        for filename in fnmatch.filter(filenames, ext):
            matches.append(os.path.join(root, filename))
    return matches

def recursive_glob2(dir_list, ext='*.xml'):

    filelist=[]
    for adir in dir_list:
        filelist.extend(recursive_glob(adir))

    return filelist

# python seismic/events.py /g/data/ha3/events_xmls_test /g/data/ha3/fxz547/travel_time_tomography/new_events20180516
if __name__ =="__main__":

    xmldir_list = sys.argv[1:]
    xmlfiles = recursive_glob2(xmldir_list)

    print ("number of events files = ", len(xmlfiles))

    print ("First file: ", xmlfiles[0])
    print ("Last file: ", xmlfiles[-1])
