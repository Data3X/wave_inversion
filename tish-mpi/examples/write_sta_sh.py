import math

data = '20240507'
nrow = 50
ncol = 50
num = str(nrow*ncol)

source_lon = 116
source_lat = 39
depth = 20
radius = 6371 - depth
lon_s = 114
lon_e = 120
lat_s = 37
lat_e = 43
mrr = 0
mrt = 1
mrp = 0
mtt = 0
mtp = 0
mpp = 0

filename = 'test2.inf'
keyword = 'c parameter for the source'
keyline = 0
with open(filename,'r') as f:
    lines = f.readlines()
    length = len(lines)
    for i in range(length):
        if keyword in lines[i]:
            keyline = i
            break
del lines[i+1:]

with open(filename,'w') as f:
    f.writelines(lines)
    f.write(' {} {} {} r0(km), lat, lon (deg)\n'.format(radius,source_lat,source_lon))
    f.write(' {} {} {} {} {} {} mt (Mrr, Mrt, Mrp, Mtt, Mtp, Mpp) (1.e25 dyne cm)\n'.format(mrr,mrt,mrp,mtt,mtp,mpp))
    f.write('c parameter for the station\n')
    f.write(' {} nr\n'.format(num))
    for i in range(nrow):
        for j in range(ncol):
            lon_sta = str(round(lon_s + j*(lon_e-lon_s)/(ncol-1),3))
            lat_sta = str(round((lat_s + i*(lat_e-lat_s)/(nrow-1)),3))
            f.write(' {} {} lat,lon (deg)\n'.format(lat_sta,lon_sta))
    f.write('c parameter for the data file\n')
    for i in range(nrow):
        for j in range(ncol):
            k = i*ncol+j+1
            name = str(k) + '.' + data + '.SH.spc'
            f.write('spcsac/example/sacfile/event2/spcfile/{}\n'.format(name))
    f.write('c\nend')