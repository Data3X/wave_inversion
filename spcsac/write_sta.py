import math

data = '20240402'
nrow = 50
ncol = 50
num = str(nrow*ncol)
lon_s = 116.20
lon_e = 116.55
lat_s = 39.75
lat_e = 40.05

with open('station.txt','w') as f:
    for i in range(nrow):
        for j in range(ncol):
            lon_sta = round(lon_s + i*(lon_e-lon_s)/(ncol-1),3)
            lat_sta = round((lat_s + j*(lat_e-lat_s)/(nrow-1)),3)
            f.write('{} {}\n'.format(lat_sta,lon_sta))