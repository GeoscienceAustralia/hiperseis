from obspy import read_inventory

def read_inspect_invxml(path2xml):

    inv = read_inventory(path2xml)
    print (inv)

if __name__ == "__main__":

    my_inventory = "/home/feizhang/Githubz/hiperseis/tests/testdata/7D_2012_2013.xml"

    read_inspect_invxml(my_inventory)