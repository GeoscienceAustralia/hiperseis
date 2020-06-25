
from collections import defaultdict

import numpy as np
import obspy

from seismic.ASDFdatabase import FederatedASDFDataSet


def main():
  src_AU_net = '../networks/network_AU_0.xml'

  au_inv = obspy.read_inventory(src_AU_net)

  ds = FederatedASDFDataSet.FederatedASDFDataSet('/g/data/ha3/Passive/SHARED_DATA/Index/asdf_files.txt')
  ds_au = ds.get_stations(obspy.UTCDateTime('1900-01-01'), obspy.UTCDateTime('2100-01-01'), network='AU')
  au_stns = defaultdict(dict)
  for row in ds_au:
    au_stns[row[1]].update({row[3]: row[4:6]})

  print("Found {} stations in asdf dataset".format(len(au_stns)))

  au_list = sorted(list(au_stns.keys()))
  missing = False
  not_close = False
  num_close = 0
  num_not_close = 0
  for au_stn in au_list:
    stn_inv = au_inv.select(station=au_stn)
    if not stn_inv:
      print("Station {} missing from inventory!".format(au_stn))
      missing = True
      continue
    else:
      assert len(stn_inv.networks) == 1
      for stn in stn_inv.networks[0]:
        for ch in stn.channels:
          if ch.code not in au_stns[au_stn]:
            continue
          expected_pos = np.array(au_stns[au_stn][ch.code])
          ch_pos = np.array([ch.longitude, ch.latitude])
          if not np.allclose(expected_pos, ch_pos, atol=0.025):
            print("{}.{} not close".format(au_stn, ch.code))
            num_not_close += 1
            not_close = True
          else:
            num_close += 1


  if not missing and not not_close:
    print("SUCCESS! ({} channel positions close)".format(num_close))
  else:
    print("FAILED! Errors found ({} channel positions not close)".format(num_not_close))
  # end if
# end func


if __name__ == '__main__':
  main()
# end if
