import sys
import netCDF4
import os
import datetime
import numpy as np
from giss import ncutil

nchunk = (15,18)    # j,i: lat,lon
IM = 360*20
JM = 180*20

YEAR = 2004

# Julian day of year for MODIS samples in YEAR (2004)
doys = (9, 17, 25, 33, 41, 49, 57, 65, 73, 81, 89, 97, 105, 113, 121, 129, 137, 145, 153, 161, 169, 177, 185, 193, 201, 209, 217, 225, 233, 241, 249, 257, 265, 273, 281, 289, 297, 305, 313, 321, 329, 337, 345, 353, 361)

# Copy these variables "plain vanilla"
# All other variables will be copied with chunking set up
vanilla_vars = {'lon', 'lat'}


def compress_dir(idir, odir):
    os.makedirs(odir, exist_ok=True)
    for doy in doys:
        sdoy = '%04d-%02d' % (YEAR, doy)
        date = datetime.datetime.strptime(sdoy, '%Y-%j')
        sdate = '%04d%02d%02d' % \
            (date.year, date.month, date.day)

        ifname = os.path.join(idir, 'Alb_soil_a.%s.006' % sdate)
        ofname = os.path.join(odir, 'Alb_soil_a_%04d_%03d.nc' % (YEAR,doy))

        compress_ncfile(ifname, ofname)
 
def compress_ncfile(ifname, ofname):

    print('---------------- Compressing to {}'.format(ofname))
    with open(ifname, 'rb') as fin:
        bytes = fin.read(JM*IM*4)
        print(len(bytes))
        val1 = np.frombuffer(bytes, dtype='f4')

    # Convert to 2D
    val2 = val1.reshape((JM,IM))

    # Flip it
    val = np.zeros(val2.shape)
    print('shape', val2.shape)
    print('JM IM', JM, IM)
    for j in range(0,JM):
        val[JM-j-1,:] = val2[j,:]
        

    with netCDF4.Dataset(ofname, 'w') as nc:
        nc.createDimension('lat', JM)
        nc.createDimension('lon', IM)

        chunk_size = (round(JM/nchunk[0]), round(IM/nchunk[1]))
        ncv = nc.createVariable('NIR', 'f', ('lat', 'lon'),
            chunksizes=chunk_size, shuffle=True, zlib=True,
            fill_value=-1.e30)
        ncv[:] = val

compress_dir('NIR_New/raw', 'NIR_New')
