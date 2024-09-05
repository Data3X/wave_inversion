import math
import numpy as np
import re

data = '200710310330A'
nrow = 46
ncol = 61
num = str(nrow*ncol)

source_lon = 116
source_lat = 39
depth = 20
radius = 6371 - depth
mrr = 38.2
mrt = 33.9
mrp = -72.1
mtt = -45.5
mtp = -10.4
mpp = 50.4

filename = 'neareq.inf'
outev = 'event3'
keyword = 'c parameter for the source'
keyline = 0
stations = [[] for i in range(ncol*nrow)]
with open('STATIONS_CS.txt','r') as station:
    lines = station.readlines()
    for i in range(ncol*nrow):
        line = lines[i].split(' ')
        stations[i].append(line[2])
        stations[i].append(line[3])

with open(filename,'r') as f:
    lines = f.readlines()
    length = len(lines)
    for i in range(length):
        if keyword in lines[i]:
            keyline = i
            break
try:    
    del lines[i+1:]
except:
    pass
        
with open(filename,'w') as f:
    f.writelines(lines)
    f.write(' {} {} {} r0(km), lat, lon (deg)\n'.format(radius,source_lat,source_lon))
    f.write(' {} {} {} {} {} {} mt (Mrr, Mrt, Mrp, Mtt, Mtp, Mpp) (1.e25 dyne cm)\n'.format(mrr,mrt,mrp,mtt,mtp,mpp))
    f.write('c parameter for the station\n')
    f.write(' {} nr\n'.format(num))
    for i in range(nrow*ncol):
        f.write('{} {} lat,lon (deg)\n'.format(stations[i][0],stations[i][1]))
    f.write('c parameter for the data file\n')
    for i in range(nrow*ncol):
        row = nrow -1 - i%nrow
        col = i//nrow
        k = row*ncol + col + 1
        name = str(k) + '.' + data + '.SH.spc'
        f.write('spcsac/example/events/' + outev + '/spcfile/{}\n'.format(name))
    f.write('c\nend')