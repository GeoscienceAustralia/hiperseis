import matplotlib.pyplot as plt
import numpy as np

loadDir='/g/data/ha3/rlt118/neural-datasets/autopicks/'

with open(loadDir+'ensemble.s.verified.txt', 'r') as f:
    vpicks = f.readlines()
    vpicks = vpicks[1:]

with open(loadDir+'ensemble.s.txt', 'r') as f:
    picks = f.readlines()


snrs = []
slopes = []
for pick in picks:
    plist = pick.split()
    snr = float(plist[17])
    slope = float(plist[20])
    snrs.append(snr)
    slopes.append(slope)
snrs = np.array(snrs)
slopes = np.array(slopes)


vsnrs = []
for vpick in vpicks:
    vplist = vpick.split()
    snr = float(vplist[17])
    vsnrs.append(snr)
vsnrs = np.array(vsnrs)

n1, bins1, a = plt.hist(snrs, bins = 20, range=(0,30), color="blue", edgecolor = "black")
n2, bins2, a = plt.hist(vsnrs, bins = 20, range=(0,30), color = "green", edgecolor = "black")
plt.xlabel("SNR")
plt.ylabel("Number")
plt.show()
bwidth = bins1[1] - bins1[0]
barcent = bins1 + bwidth/2
barcent = barcent[0:-1]
plt.bar(barcent,n2/n1,width=bwidth, color = "blue", edgecolor = "black")
plt.xlabel("SNR")
plt.ylabel("Acceptance ratio")
plt.show()
# plt.scatter(snrs,slopes)
# plt.show()
