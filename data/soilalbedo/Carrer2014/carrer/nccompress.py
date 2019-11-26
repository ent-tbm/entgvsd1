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
doys = (1, 9, 17, 25, 33, 41, 49, 57, 65, 73, 81, 89, 97, 105, 113, 121, 129, 137, 145, 153, 161, 169, 177, 185, 193, 201, 209, 217, 225, 233, 241, 249, 257, 265, 273, 281, 289, 297, 305, 313, 321, 329, 337, 345, 353, 361)

# Copy these variables "plain vanilla"
# All other variables will be copied with chunking set up
vanilla_vars = {'lon', 'lat'}

# Cutoff points for spectral bands
bands = (300,770,5000)
band_labels = ('VIS', 'NIR')

def get_ifiles(idir):
    ifiles = []
    for doy in doys:
        sdoy = '%04d-%02d' % (YEAR, doy)
        date = datetime.datetime.strptime(sdoy, '%Y-%j').date()
        sdate = '%04d%02d%02d' % \
            (date.year, date.month, date.day)

        ifname = os.path.join(idir, 'Alb_soil_a.%s.006' % sdate)
        ifiles.append((date, ifname))


    return ifiles


FillValue = -1.e30

def create_ncfile(ofname, bands, dates):
    with netCDF4.Dataset(ofname, 'w') as ncout:
        ncout.createDimension('lat', JM)
        nclat = ncout.createVariable('lat', 'd', ('lat',))
        nclat.long_name = 'latitude'
        nclat.units = 'degrees north'
        dy = 180. / JM
        nclat[:] = np.array([(i+.5)*dy-90. for i in range(0,JM)])

        ncout.createDimension('lon', IM)
        nclon = ncout.createVariable('lon', 'd', ('lon',))
        nclon.long_name = 'longitude'
        nclon.units = 'degrees east'
        dx = 360. / IM
        nclon[:] = np.array([(i+.5)*dx-180. for i in range(0,IM)])

        ncout.createDimension('dates', len(dates))
        ncout.createDimension('datelen', 10)    # Length of date string YYYY-MM-dd
        nbands = len(bands)-1
        ncout.createDimension('bands', nbands)
        ncout.createDimension('bands_plus1', nbands+1)

        ncout.title = 'Albedo of Soil'
        ncout.history = 'Aug 2019: Compiled from non-NetCDF files provided by Carrer'
        ncout.creator_name = 'Elizabeth Fischer, Dominique Carrer'
        ncout.creator_email = 'elizabeth.fischer@columbia.edu, dominique.carrer@meteo.fr'
        ncout.geospatial_lat_min = -90
        ncout.geospatial_lat_max = 90
        ncout.geospatial_lon_min = -180
        ncout.geospatial_lon_max = 180


        # Store bands
#        ncv = ncout.createVariable('banddiv', 'd', ('bands_plus1',))
#        ncv.long_name = 'Definition of spectral bands'
#        ncv.units = 'nm'
#        ncv[:] = bands

#        ncout.createDimension('bands.strlen', 18)
#        ncv = ncout.createVariable('bands', 'c', ('bands','bands.strlen'))
#        ncv.long_name = 'Spectral range of each band'
#        ncv.units = 'nm'
        sbands = []
        band_ranges = []
        for i in range(0,nbands):
            sbands.append('%s %04d - %04d nm' % (band_labels[i],bands[i],bands[i+1]))
            band_ranges.append((float(bands[i]),float(bands[i+1])))

#        ncv[:] = netCDF4.stringtochar(np.array(sbands, dtype='S'))


#        ncout.createDimension('band_labels.strlen', 3)
#        ncv = ncout.createVariable('band_labels', 'c', ('bands', 'band_labels.strlen'))
#        ncv[:] = netCDF4.stringtochar(np.array(band_labels, dtype='S'))

        # Store dates as integers
        dtbase = datetime.date(1970,1,1)
        ncv = ncout.createVariable('dates', 'i', ('dates',))
        ncv.setncattr('units', 'days since 1970-01-01')
        ncv.setncattr('description', 'Date of each observation')
        ncv[:] = [(dt - dtbase).days for dt in dates]

        # Store dates as strings
        ncv = ncout.createVariable('sdates', 'c', ('dates','datelen'))
        ncv.setncattr('description', 'Date of each observation, as string')
        ncv[:] = netCDF4.stringtochar(np.array([dt for dt in dates], dtype='S'))

        chunk_size = (1, round(JM/nchunk[0]), round(IM/nchunk[1]))
        for blab,sband,band_range in zip(band_labels,sbands,band_ranges):
            ncv = ncout.createVariable('soilalb_{}'.format(blab), 'f',
                ('dates', 'lat', 'lon'),
                chunksizes=chunk_size, shuffle=True, zlib=True,
                fill_value=FillValue)
            ncv.band_description = sband
            ncv.band_range = band_range
            ncv.units = '1'



def compress_ncfiles(ifiles, iband, ofname):

    # Check that all ifiles exist
    nmissing = 0
    for _,ifile in ifiles:
        if (not os.path.exists(ifile)):
            sys.stderr.write('Missing file: {}\n'.format(ifile))
            nmissing += 1
    if nmissing > 0:
        sys.exit(-1)

    with netCDF4.Dataset(ofname, 'a') as ncout:

        ncv = ncout.variables['soilalb_{}'.format(band_labels[iband])]

        for idt,(dt,ifname) in enumerate(ifiles):
            print('Reading {}'.format(ifname))
            with open(ifname, 'rb') as fin:
                bytes = fin.read(JM*IM*4)
                val1 = np.frombuffer(bytes, dtype='f4')

            # Convert to 2D
            val2 = val1.reshape((JM,IM))

            # Flip it
            val = np.zeros(val2.shape)
            for j in range(0,JM):
                val[JM-j-1,:] = val2[j,:]

            # Convert NetCDF FillValue
            val[val == 999] = FillValue

            # Store it
            ncv[idt,:,:] = val


YEAR = 2004
dates = [
    datetime.datetime.strptime('%04d-%02d' % (YEAR, doy), '%Y-%j').date()
    for doy in doys]

create_ncfile('carrer.nc', bands, dates)

compress_ncfiles(get_ifiles('VIS_2016'), 0, 'carrer.nc')
compress_ncfiles(get_ifiles('NIR_New'), 1, 'carrer.nc')

