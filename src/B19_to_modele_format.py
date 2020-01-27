# Converts ModelE input data file to form usable by ModelE
# This fixes two things:
#   1. Reworks the data to the NetCDF scheme expected by ModelE
#   2. Writes NetCDF3 output files
#
# Author: Elizabeth Fischer
#
import netCDF4
import re
import os
import itertools
import collections

def decode_strs(ncin, vname):
    """Reads a bunch of strings and returns as a list"""
    return [x.decode().strip() for x in netCDF4.chartostring(ncin.variables[vname][:])]

imonths = dict((m,im+1) for im,m in
    enumerate("Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec".split('|')))
smonths = dict((im+1,m) for im,m in
    enumerate("Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec".split('|')))



File = collections.namedtuple('File', ('ftype', 'name', 'imonth', 'mname'))
fileRE = re.compile(r'(.*)_((lai)|(laimax)|(hgt)|(lc))_(.*)\.nc')
monthRE = re.compile(r'(.*)_(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)_(.*)\.nc')
def snoop_A10_dir(idir):
    """Snoops output of B18_modele to determine the different relevant files in it"""
    files = []
    for fname in os.listdir(idir):
        if fname.endswith('.nc3'):
            continue
        match = fileRE.match(fname)
        if match is not None:
            ftype = match.group(2)
            name = os.path.join(idir, '{}_{}_{}'.format(match.group(1),match.group(2),match.group(7)))
            match = monthRE.match(fname)
            if match is None:
                imonth = None
                mname = None
            else:
                imonth = imonths[match.group(2)]
                mname = os.path.join(idir, '{}_monthly_{}'.format(match.group(1),match.group(3)))
            files.append(File(ftype,name,imonth,mname))

    annual_files = list(file for file in files if file.imonth is None)
    month_groups = list((k,list(g)) for k,g in itertools.groupby(sorted(
        (file for file in files if file.imonth is not None), key=lambda f: (f.mname,f.imonth)),
        key=lambda f: f.mname))

    return annual_files, month_groups
# ----------------------------------------------
def copydef_base(ncin, ncout):
    # Copy lat and lon vars
    for vname in ('lat', 'lon'):
        ncout.createDimension(vname, len(ncin.dimensions[vname]))
        ncv = ncout.createVariable(vname, 'f', (vname,))
        ncv.units = ncin.variables[vname].units

    # Copy global attributs
    for attr in ncin.ncattrs():
        if attr != 'long_name':
            ncout.setncattr(attr, ncin.getncattr(attr))

def copy_base(ncin, ncout):
    for vname in ('lat', 'lon'):
        ncout.variables[vname][:] = ncin.variables[vname][:]
# ----------------------------------------------
def convert_annual(ivname, ifname, ofname):
    print('------------------- {}'.format(ifname))
    with netCDF4.Dataset(ifname) as ncin:
        with netCDF4.Dataset(ofname, 'w', format="NETCDF3_CLASSIC") as ncout:
            # ------------- Define
            copydef_base(ncin, ncout)
            iv = ncin.variables[ivname]

            ncout.long_name = iv.long_name
            units = iv.units
            _FillValue = iv._FillValue

            layers = decode_strs(ncin, 'lctype')
            long_names = decode_strs(ncin, 'lctype_longname')

            for i,(layer,long_name) in enumerate(zip(layers,long_names)):
                ncv = ncout.createVariable(layer, 'f', ('lat', 'lon'), fill_value=_FillValue)
                ncv.units = units
                ncv.long_name = '{} - {}'.format(i+1, long_name)

            # -------------- Copy
            for i,layer in enumerate(layers):
                ncout.variables[layer][:] = ncin.variables[ivname][i,:]

# ----------------------------------------------
def convert_annual(ivname, ifname, ofname, prefix=''):
    print('------------------- {}'.format(ifname))
    with netCDF4.Dataset(ifname) as ncin:
        with netCDF4.Dataset(ofname, 'w', format="NETCDF3_CLASSIC") as ncout:
            # ------------- Define
            copydef_base(ncin, ncout)
            iv = ncin.variables[ivname]
            print('variable ', ivname)
            ncout.long_name = iv.long_name
            units = iv.units
            _FillValue = iv._FillValue

            layers = decode_strs(ncin, 'lctype')
            long_names = decode_strs(ncin, 'lctype_longname')

            for i,(layer,long_name) in enumerate(zip(layers,long_names)):
                ncv = ncout.createVariable(prefix+layer, 'f', ('lat', 'lon'), fill_value=_FillValue)
                ncv.units = units
                ncv.long_name = '{} - {}'.format(i+1, long_name)

            # -------------- Copy
            copy_base(ncin, ncout)
            for i,layer in enumerate(layers):
                ncout.variables[prefix+layer][:] = ncin.variables[ivname][i,:]

# ----------------------------------------------
def convert_monthly(ivname, mname, files, ifiles, ofiles):
    """
    ifiles: OUT
        Append files we opened for input to this list
    ofoiles: OUT
        Append files we opened for output to this list"""

    print('------------------- {}'.format(mname))
    with netCDF4.Dataset(mname + '.nc3', 'w', format="NETCDF3_CLASSIC") as ncout:
        ofiles.append(mname + '.nc3')

        # ----------- Define
        with netCDF4.Dataset(files[0].name+'.nc') as ncin:
            copydef_base(ncin, ncout)
            ncout.createDimension('time', len(files))
            nctime = ncout.createVariable('time', 'i', ('time',))
            nctime.units = "Month Sequence of Climatology"

            iv = ncin.variables[ivname]
            match = re.match(r'(.*), {}'.format(smonths[files[0].imonth]), iv.long_name)
            if match is None:
                ncout.long_name = iv.long_name
            else:
                ncout.long_name = match.group(1) + ', Monthly'
            units = iv.units
            _FillValue = iv._FillValue

            layers = decode_strs(ncin, 'lctype')
            long_names = decode_strs(ncin, 'lctype_longname')

            for i,(layer,long_name) in enumerate(zip(layers,long_names)):
                ncv = ncout.createVariable(layer, 'f', ('time', 'lat', 'lon'), fill_value=_FillValue)
                ncv.units = units
                ncv.long_name = '{} - {}'.format(i+1, long_name)

        # ------------------ Copy
        with netCDF4.Dataset(files[0].name+'.nc') as ncin:
            copy_base(ncin, ncout)

        nctime[:] = [file.imonth for file in files]
        for itime,file in enumerate(files):
            with netCDF4.Dataset(file.name+'.nc') as ncin:
                ifiles.append(file.name+'.nc')
                for ilayer,(layer,long_name) in enumerate(zip(layers,long_names)):
                    ncout.variables[layer][itime,:] = ncin.variables[ivname][ilayer,:]

# ----------------------------------------------

convert_fn = {
    'lc' : lambda ifname,ofname: convert_annual('lc', ifname, ofname),
    'laimax' : lambda ifname,ofname: convert_annual('laimax', ifname, ofname),
    'lai' : lambda ifname,ofname: convert_monthly('lai', ifname, ofname),
    'hgt' : lambda ifname,ofname: convert_annual('hgt', ifname, ofname, prefix='hgt_')
}

OUTPUTS_DIR = '../outputs'

def main():
    annual_files, month_group = snoop_A10_dir(os.path.join(OUTPUTS_DIR, 'modelE'))

    ifiles = list()
    ofiles = list()
    for file in annual_files:
        if file.ftype in convert_fn:
            ifname = file.name+'.nc'
            ofname = file.name+'.nc3'
            convert_fn[file.ftype](ifname, ofname)

    for mname,files in month_group:
        convert_monthly('lai', mname,files, ifiles,ofiles)

    # Write the .mk file, to be compatible with the rest of this suite
    exe = 'B19_to_modele_format'
    with open(os.path.join(OUTPUTS_DIR, exe+'.mk'), 'w') as out:
        out.write('{}_INPUTS = \\\n'.format(exe))
        for fname in ifiles:
            out.write('    {} \\\n'.format(fname))
        out.write('\n')

        out.write('{}_OUTPUTS = \\\n'.format(exe))
        for fname in ofiles:
            out.write('    {} \\\n'.format(fname))
        out.write('\n')

main()
