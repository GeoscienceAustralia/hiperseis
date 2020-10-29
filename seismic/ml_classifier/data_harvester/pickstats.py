import csv, bisect

with open('/g/data/ha3/rlt118/neural-datasets/categoriser-teleseismic/smallset-2/picklog.csv') as f:
    picklog = csv.reader(f)
    srted = []
    for pick in picklog:
        dist = float(pick[-1].strip())
        bisect.insort(srted,dist)

with open('/g/data/ha3/rlt118/neural-datasets/categoriser-teleseismic/smallset-2/sortdist.txt','w') as g:
    for flt in srted:
        g.write(str(flt)+'\n')


for idx, flt in enumerate(srted[::-1]):
    if flt <= 80:
        print(idx)
        break

