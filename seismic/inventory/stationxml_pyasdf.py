"""
GA AusArray has extra metadata (json/csv-formatted strings) in stationXML file.
Ingest into ASDF h5, extract the stationxml and write out to a file.
Compare the original statoinxml with the extracted output xml to ensure they are identical.

Fei.zhang@ga.gov.au
2020-07-30
"""
import obspy
import pyasdf

if __name__ == "__main__":

    print("pyasdf.__version__", pyasdf.__version__, pyasdf.__file__)
    print("obspy.__version__", obspy.__version__, obspy.__file__)

    filename = "./OA.CF28.xml"  # input xml file with extra metadata in it already

    # Read directly.
    inv_original = obspy.read_inventory(filename)

    # Write to ASDF file, and read again.
    asdf_file = "test.h5"
    with pyasdf.ASDFDataSet(asdf_file) as ds:
        ds.add_stationxml(filename)

    with pyasdf.ASDFDataSet(asdf_file) as ds:
        inv_new = ds.waveforms["OA.CF28"].StationXML

        inv_new.write('OA.CF28_new.xml', format='STATIONXML', nsmap={
                      'GeoscienceAustralia': 'https://github.com/GeoscienceAustralia/hiperseis/xmlns/1.0'})

    assert inv_original == inv_new

    assert (
        inv_original.networks[0].stations[0].extra
        == inv_new.networks[0].stations[0].extra
    )

    print(inv_original.networks[0].stations[0].extra)
    print(inv_new.networks[0].stations[0].extra)

    # Now, read the file 'OA.CF28_new.xml' try to get the extra metadata. Should be OK.
    inv_new2 = obspy.read_inventory('OA.CF28_new.xml')
    print(inv_new2.networks[0].stations[0].extra)

    # Using diff to compare the origial input xml file "OA.CF28.xml" to the output 'OA.CF28_new.xml'.
    # They should be essentially identical, except some minor diffs.
