import sys
import netCDF4
import os
from giss import ncutil

nchunk = (15,18)

# Copy these variables "plain vanilla"
# All other variables will be copied with chunking set up
vanilla_vars = {'lon', 'lat'}

def compress_ncfile(fin, fout):

    print('---------------- Compressing to {}'.format(fout))
    with netCDF4.Dataset(fin, 'r') as ncin:
        with netCDF4.Dataset(fout, 'w') as ncout:

            im = len(ncin.dimensions['lon'])
            jm = len(ncin.dimensions['lat'])

            ncc = ncutil.copy_nc(ncin, ncout)

            compress_vars = [x for x in ncin.variables if x not in vanilla_vars]

            chunk_size = (round(jm/nchunk[0]), round(im/nchunk[1]))
            ncc.define_vars(vanilla_vars)
            ncc.define_vars(compress_vars, chunksizes=chunk_size, shuffle=True, zlib=True, fill_value=-1.e30)

#            ncout.define_vars(chunksizes=shuffle=True)

#            ncc.copy_data()


            for vname in vanilla_vars:
                ncc.copy_var(vname, vname)

            for vname in compress_vars:
                print('Compressing {}'.format(vname))
                ivar = ncin.variables[vname]
                ovar = ncout.variables[vname]
                for jchunk in range(0,nchunk[0]):
                    print('    ', end='')
                    for ichunk in range(0,nchunk[1]):
                        print('.', end='')
#                        print('    chunk ({}, {})'.format(jchunk,ichunk))
                        sys.stdout.flush()
                        ji0 = (jchunk*chunk_size[0], ichunk*chunk_size[1])
#                        print(ji0[0],ji0[0]+chunk_size[0], ji0[1],ji0[1]+chunk_size[1])

                        val = ivar[ji0[0]:ji0[0]+chunk_size[0], ji0[1]:ji0[1]+chunk_size[1]]
                        ovar[ji0[0]:ji0[0]+chunk_size[0], ji0[1]:ji0[1]+chunk_size[1]] = val
                    print('')




#compress_ncfile('SW_Alb_soil_yearly.006.1kmx1km.nc', 'x.nc')

for fname in os.listdir('nc3'):
    compress_ncfile(os.path.join('nc3',fname), fname)
