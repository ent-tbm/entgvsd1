import netCDF4

exclude = {'source'}

for fname0,fname1 in (
    ('EntGVSD_v1.1_Monfreda_crops_5min.nc', 'EntGVSD_v1.1_Monfreda_crops_1km.nc'),
    ('EntGVSD_v1.1_Monfreda_crops_5min_norm.nc', 'EntGVSD_v1.1_Monfreda_crops_1km_norm.nc'),
    ('EntGVSD_v1.1_Monfreda_crops_5min.nc', 'EntGVSD_v1.1_Monfreda_crops_1km_forplot.nc'),
    ('EntGVSD_v1.1_Monfreda_crops_5min_norm.nc', 'EntGVSD_v1.1_Monfreda_crops_1km_norm_forplot.nc')):


    print('------------ {} {}'.format(fname0, fname1))
    with netCDF4.Dataset(fname0) as nc0:
        with netCDF4.Dataset(fname1, 'a') as nc1:
            for v0 in nc0.variables.values():
                v1 = nc1.variables[v0.name]
                for attr_name in v0.ncattrs():
                    if attr_name not in exclude:
                        attr_val = v0.getncattr(attr_name)
                        print('    {}.{} = {}'.format(v1.name, attr_name, attr_val))
                        v1.setncattr(attr_name, attr_val)
