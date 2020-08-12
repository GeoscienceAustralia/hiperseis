"""
GA AusArray has extra metadata (json/csv-formatted strings) in stationXML file.
Using pyasdf, we ingest the stationxml into ASDF h5, then extract the stationxml and write out to a file.
Compare the original statoinxml with the extracted output xml
and check if they are essentially identical?

There were problems before pyasdf version 0.7.1.
Should install the latest version 0.7.3 with python3.7+

Ref:
https://github.com/SeismicData/pyasdf/issues/61
https://github.com/SeismicData/pyasdf/issues/63

Fei.zhang@ga.gov.au
2020-07-30
"""
import obspy
import pyasdf

if __name__ == "__main__":

    print("pyasdf.__version__", pyasdf.__version__, pyasdf.__file__)
    print("obspy.__version__", obspy.__version__, obspy.__file__)

    input_file_2_asdf = "./OA.CF28_input2_pyasdf.xml"  # input xml file with extra metadata in it already
    output_from_asdf =  "./OA.CF28_out_pyasdf.xml"
    # Read directly.
    inv_original = obspy.read_inventory(input_file_2_asdf)

    # Write to ASDF file, and read again.
    asdf_file = "test.h5"
    with pyasdf.ASDFDataSet(asdf_file) as ds:
        ds.add_stationxml(input_file_2_asdf)

    with pyasdf.ASDFDataSet(asdf_file) as ds:
        inv_new = ds.waveforms["OA.CF28"].StationXML

        inv_new.write(output_from_asdf, format='STATIONXML', nsmap={
                      'GeoscienceAustralia': 'https://github.com/GeoscienceAustralia/hiperseis/xmlns/1.0'})

    assert inv_original == inv_new

    assert (
        inv_original.networks[0].stations[0].extra
        == inv_new.networks[0].stations[0].extra
    )

    print(inv_original.networks[0].stations[0].extra)
    print(inv_new.networks[0].stations[0].extra)

    # Now, read the xml extracted out from the ASDF file and try to access the extra metadata.
    # Should be working if the xml format is standard obspy/FSDN .
    inv_new2 = obspy.read_inventory(output_from_asdf)
    print(inv_new2.networks[0].stations[0].extra)

    # Using diff to compare the origial input xml file "OA.CF28.xml" to the output 'OA.CF28_new.xml'.
    # They should be essentially identical, except some minor diffs.
