from obspy import read_inventory
from obspy import Inventory
from obspy.core.inventory import Network
from obspy.core.util import AttribDict

def get_csv_correction_data(path_csvfile):
    """
    Read in the csv data from an input file, get the network_code, station_code, csv_data

    Args:
        path_csvfile: input csv file with data like
        net,sta,date,clock_correction
        7D,DE43,2012-11-27,1.0398489013215846
        7D,DE43,2012-11-28,0.9408504322549281
        7D,DE43,2012-11-29,0.8418519631882714
        7D,DE43,2012-11-30,0.7428534941216148
        7D,DE43,2012-12-01,0.6438550250549583

    Returns: (network_code, station_code, csv_data)

    """
    with open(path_csvfile,"r") as csvfid:
        all_lines = csvfid.readlines()

    print("The csv file length", len(all_lines))
    print("the first line", all_lines[0])
    print("the last  line", all_lines[-1])

    line2 = all_lines[1]
    print("line2=",line2)

    items =line2.split(",")

    network_code = items[0]
    station_code = items[1]

# fxz547@vdi-n25 /g/data/ha3/Passive/SHARED_DATA/GPS_Clock/corrections
# $ head 7D.DE43_clock_correction.csv
#     csv_data = """
#     net,sta,date,clock_correction
#     7D,DE43,2012-11-27,1.0398489013215846
#     7D,DE43,2012-11-28,0.9408504322549281
#     7D,DE43,2012-11-29,0.8418519631882714
#     7D,DE43,2012-11-30,0.7428534941216148
#     7D,DE43,2012-12-01,0.6438550250549583
#     7D,DE43,2012-12-02,0.5448565559883017
#     7D,DE43,2012-12-03,0.445858086921645
#     7D,DE43,2012-12-04,0.3468596178549885
#
#     """
#
#     network_code = "7D"
#     station_code = "DE43"

    return (network_code,station_code, all_lines)

def add_gpscorrection_into_stationxml(csv_file, path2_myxml):

    ns = "https://github.com/GeoscienceAustralia/hiperseis/xmlns/1.0"

    (net,sta,csv_data) = get_csv_correction_data(csv_file)

    # path2_myxml = "/home/feizhang/Githubz/hiperseis/tests/testdata/7D_2012_2013.xml"
    my_inv = read_inventory(path2_myxml)

    # https://docs.obspy.org/packages/autogen/obspy.core.inventory.inventory.Inventory.select.html#obspy.core.inventory.inventory.Inventory.select

    selected_inv = my_inv.select(station=sta)

    print(selected_inv)


    my_tag = AttribDict()
    my_tag.namespace = ns
    my_tag.value = csv_data


    selected_inv.networks[0].stations[0].extra = AttribDict()
    selected_inv.networks[0].stations[0].extra.gpsclockcorrection =my_tag

    stationxml_with_csv ='modified_inventory_%s.xml'%sta
    selected_inv.write(stationxml_with_csv, format='STATIONXML',
               nsmap={'GeoscienceAustralia': 'https://github.com/GeoscienceAustralia/hiperseis/xmlns/1.0'})

    #my_inv.write('modified_inventory.xml', format='STATIONXML')

    return stationxml_with_csv


def extract_csvdata(path2xml):
    """
    Read and extract the csv data from an inventory file

    Args:
        path2xml: path_to_stationxml

    Returns:
        csv_str

    """

    new_inv = read_inventory(path2xml)

    csv_str = new_inv.networks[0].stations[0].extra.gpsclockcorrection.value

    return csv_str


if __name__ == "__main__":

    #VDI   time_correction_csvfile = "/g/data/ha3/Passive/SHARED_DATA/GPS_Clock/corrections/input_file"
    time_correction_csvfile ="/home/feizhang/Githubz/hiperseis/tests/testdata/corrections/7D.DE43_clock_correction.csv"
    time_correction_csvfile = "/home/feizhang/Githubz/hiperseis/tests/testdata/corrections/7D.CZ40_clock_correction.csv"

    my_inventory = "/home/feizhang/Githubz/hiperseis/tests/testdata/7D_2012_2013.xml"

    out_xml = add_gpscorrection_into_stationxml(time_correction_csvfile, my_inventory )

    csvstr= extract_csvdata(out_xml)

    print(csvstr)
    print(type(csvstr))