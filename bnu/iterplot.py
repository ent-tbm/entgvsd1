import netCDF4
import giss.basemap
import giss.plot
import matplotlib.pyplot
import mpl_toolkits.basemap
import sys
import numpy as np
import os
import collections
import psutil
import fcntl
import collections

PlotSpec = collections.namedtuple('PlotSpec', ('ifname', 'vname', 'layer', 'algo', 'ofname'))
Region = collections.namedtuple('Region',
    ('name',        # Descriptive name of the region
    'center_ll'))   # (lat, lon) of the center of the region

# Convert (deg,min,sec) into decimal lat/lon
def _N(deg,min,sec):
    return deg + min/60. + sec/3600.
def _S(deg,min,sec):
    return -(deg + min/60. + sec/3600.)
def _E(deg,min,sec):
    return deg + min/60. + sec/3600.
def _W(deg,min,sec):
    return - (deg + min/60. + sec/3600.)

regions = {
    'ks' : Region('Kansas', (39., -98.)),    # (lat,lon)
    'afk' : Region('Sub-Saharan Africa', (_N(11,16,00.9), _E(18,59,44.2))),
    'eur' : Region('Europe crop belt', (_N(50,35,58.9), _E(19,50,00.1))),
    'bra' : Region('Brazil crops/shrubland', (_S(17,25,34.2), _W(53,0,40.5))),
}


def plot_global(ps):
    """region.center_ll:
        (lat, lon) of center of region to plot"""
    print('Plotting {}'.format(ps.ifname))

    with netCDF4.Dataset(ps.ifname) as nc:

        nlon = nc.variables['lon'].shape[0]
        nlat = nc.variables['lat'].shape[0]

        screenres = (1800,2880)    # y,x
        dotpitch = 1./220    # 220 dpi monitor

        # Read lons and lats (cell centers)
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]

        # Get data variable name
        ncvar = nc.variables[ps.vname]
        if ps.layer < 0:
            val = ncvar[:]
        else:
            val = ncvar[ps.layer,:]
        val[val==0] = np.nan
        plotter = giss.plot.LonLatPlotter(lon, lat, boundaries=False)
     
        basemap = giss.basemap.global_map()

        # One plot per page (figsize must be in inches)
        figure = matplotlib.pyplot.figure(figsize=(screenres[1]*dotpitch, screenres[0]*dotpitch))

        ax = figure.add_subplot(111)

        cb_args=dict()
        if ps.ifname.endswith('_lc_lr.nc') or ps.ifname.endswith('brightratio_lr.nc'):
            plot_args=dict(vmin=0.,vmax=1.0)

        elif ps.ifname.endswith('_lai_lr.nc'):
            plot_args=dict(vmin=0.,vmax=7.0)
            plot_args['cmap'] = giss.plot.read_cpt('NYT_drought14b.cpt').cmap
        else:
            raise ValueError('Unknown plot type (lc vs lai)')

        _title = nc.title if hasattr(nc,'title') else os.path.split(ps.ifname)[1]
        plt = giss.plot.plot_var(ax=ax, basemap=basemap,
            show=False,
            plotter=plotter,
            val=val,
            title='{}\n{}'.format(_title, ncvar.long_name),
            plot_args=plot_args, cb_args=cb_args)

        basemap.drawcoastlines(linewidth=.1)

        # draw parallels and meridians.
        # label parallels on right and top
        # meridians on bottom and left
        # labels = [left,right,top,bottom]
        meridians = np.arange(-180.,180.,30.)
        basemap.drawmeridians(meridians,labels=[False,False,False,True])

        parallels = np.arange(-90.,90.,30.)
        basemap.drawparallels(parallels,labels=[True,False,False,False])

        # Save to a file as png
        figure.tight_layout(pad=2.0)  # https://stackoverflow.com/questions/4042192/reduce-left-and-right-margins-in-matplotlib-plot
        print('    Writing {}'.format(ps.ofname))
        odir = os.path.split(ps.ofname)[0]
        os.makedirs(odir, exist_ok=True)
        tmp_fname = os.path.join(odir, '_tmp.png')
        figure.savefig(tmp_fname, dpi=1./dotpitch, transparent=False)
        os.rename(tmp_fname, ps.ofname)   # Write atomically

def plot_region(ps):
#ifname, ofname, region):
    """region.center_ll:
        (lat, lon) of center of region to plot"""
    region = regions[ps.algo]

    if ps.ifname.endswith('_err.nc'):
        return

    print('Plotting {} ({})'.format(ps.ifname, region.name))

    with netCDF4.Dataset(ps.ifname) as nc:

        nlon = nc.variables['lon'].shape[0]
        nlat = nc.variables['lat'].shape[0]

        screenres = (1800,2880)    # y,x
        dotpitch = 1./220    # 220 dpi monitor

        # Convert region center from geographic lon/lat to i/j
#        region.center_ll = (39., -98)   # (lat,lon)
    #    region.center_ll = (-180+dxy[0]*1.5001,90.)   # (lon, lat)
        center_ji = (    # Zero-based indexing (lat,lon)
            int(round((90.+region.center_ll[0]) * (nlat-1) *  (1./180.))),
            int(round((180.+region.center_ll[1]) * (nlon-1) * (1./360.))))

        # Determine region bounds
        # region_size = (int(screenres[0]*.8),int(screenres[1]*.95))    # Size of region, in gridcells (y,x)
        region_size = (1320,2600)  # Size of region to plot, in gridcells (y,x)
        _base = (center_ji[0]-region_size[0]//2, center_ji[1]-region_size[1]//2)
        _top = (_base[0]+region_size[0]), _base[1]+region_size[1]
        region_ji = (_base, _top)

        if (_base[0] < 0 or _base[1] < 0 or _top[0] > nlat or _top[1] > nlon):
            raise ValueError('Coordinates out of bounds', _base, _top, region.center_ll)

        # Read lons and lats (cell centers)
        lon = nc.variables['lon'][_base[1]:_top[1]]
        lat = nc.variables['lat'][_base[0]:_top[0]]

        # Convert to boundaries
        lonb = np.zeros(len(lon)+1)
        hdlon = .5 * 360. / nlon
        lonb[0:-1] = lon - hdlon
        lonb[-1] = lon[-1] + hdlon

        latb = np.zeros(len(lat)+1)
        hdlat = .5 * 180. / nlat
        latb[0:-1] = lat - hdlat
        latb[-1] = lat[-1] + hdlat


        # Get data variable name
        ncvar = nc.variables[ps.vname]
        if ps.layer < 0:
            val = ncvar[_base[0]:_top[0],  _base[1]:_top[1]]
        else:
            val = ncvar[ps.layer, _base[0]:_top[0],  _base[1]:_top[1]]
        val[val==0] = np.nan
        plotter = giss.plot.LonLatPlotter(lonb, latb, boundaries=True)
     
        print('    Region ({:.2f}, {:.2f}) -- ({:.2f}, {:.2f})'.format(lonb[0],latb[0],lonb[-1],latb[-1]))


        # Use a custom basemap
        #resolution='l',
        basemap = mpl_toolkits.basemap.Basemap(
            projection='laea',\
            lat_0=region.center_ll[0], lon_0=region.center_ll[1],
            llcrnrlon=lonb[0], llcrnrlat=latb[0],
            urcrnrlon=lonb[-1], urcrnrlat=latb[-1])

        # One plot per page (figsize must be in inches)
        figure = matplotlib.pyplot.figure(figsize=(screenres[1]*dotpitch, screenres[0]*dotpitch))

        ax = figure.add_subplot(111)

        cb_args=dict()
        if ps.ifname.endswith('_lc.nc') or ps.ifname.endswith('brightratio.nc'):
            plot_args=dict(vmin=0.,vmax=1.0)
        elif ps.ifname.endswith('_lai.nc'):
            plot_args=dict(vmin=0.,vmax=7.0)
            plot_args['cmap'] = giss.plot.read_cpt('NYT_drought14b.cpt').cmap
        else:
            raise ValueError('Unknown plot type (lc vs lai)')

        _title = nc.title if hasattr(nc,'title') else os.path.split(ps.ifname)[1]
        plt = giss.plot.plot_var(ax=ax, basemap=basemap,
            show=False,
            plotter=plotter,
            val=val,
            title='{}\n{}: ({})'.format(_title, ncvar.long_name, region.name),
            plot_args=plot_args, cb_args=cb_args)

        # draw parallels and meridians.
        # label parallels on right and top
        # meridians on bottom and left
        # labels = [left,right,top,bottom]
        meridians = np.arange(int(lonb[0])-1., int(lonb[-1])+1., 4.)
        basemap.drawmeridians(meridians,labels=[False,False,False,True])

        parallels = np.arange(int(latb[0])-1., int(latb[-1])+1., 4.)
        basemap.drawparallels(parallels,labels=[True,False,False,False])

        # Save to a file as png
        figure.tight_layout(pad=2.0)  # https://stackoverflow.com/questions/4042192/reduce-left-and-right-margins-in-matplotlib-plot
        print('    Writing {}'.format(ps.ofname))
        odir = os.path.split(ps.ofname)[0]
        os.makedirs(odir, exist_ok=True)
        tmp_fname = os.path.join(odir, '_tmp.png')
        figure.savefig(tmp_fname, dpi=1./dotpitch, transparent=False)
        os.rename(tmp_fname, ps.ofname)   # Write atomically



def iter_plot():
    for path,dirnames,filenames in os.walk('lc_lai_ent'):
        for fname in sorted(filenames):
            if fname.endswith('.nc') and not fname.endswith('_lr.nc'):
                # Determine file type
                ifname = os.path.join(path,fname)
                odir = os.path.splitext(ifname)[0]

#                if (not odir.endswith('.nc')):
#                    continue
#                yield odir
#                continue

                with netCDF4.Dataset(ifname) as nc:
                    # Determine how may plots in each variable
                    if 'layers' in nc.dimensions:
                        nlayer = len(nc.dimensions['layers'])
                        layer_names = []
                        ln = nc.variables['layer_names'][:]
                        for i in range(0,ln.shape[0]):
                            name_b = ln[i,:].tobytes().decode('utf-8').strip()
                            layer_names.append(name_b)
                    else:
                        nlayer = 0


                    # Loop through variables
                    for vname in nc.variables:
                        if vname in {'lon', 'lat', 'layer_indices', 'layer_names'}:
                            continue

                        if nlayer == 0:
                            # Plot this variable a single time
                            yield PlotSpec(odir + '_lr.nc', vname, -1, 'global',
                                os.path.join(odir, 'global.png'))
                            for rname in regions.keys():
                                yield PlotSpec(ifname, vname, -1, rname,
                                    os.path.join(odir, rname + '.png'))
                        else:
                            # Plot each layer in this variable
                            for ix in range(0,nlayer):
                                yield PlotSpec(odir + '_lr.nc', vname, ix, 'global',
                                    os.path.join(odir, layer_names[ix]+'_global.png'))
                                for rname in regions.keys():
                                    yield PlotSpec(ifname, vname, ix, rname,
                                        os.path.join(odir, '{}_{}.png'.format(layer_names[ix], rname)))

def iter_plot_todo():
    for ps in iter_plot():
        # Check to see if work is complete.
        if os.path.exists(ps.ofname) and (os.path.getmtime(ps.ifname) <= os.path.getmtime(ps.ofname)):
            continue


        # Check for ourselves if we are runnig low on memory.
        # If it's this low, probably another plotregion is running.
        # Exit with no error, so our caller doesn't keep calling us.
        avail = psutil.virtual_memory().available
        print('********** avail', avail)
        if avail < 1000000000:
            sys.exit(1)

        print('mem1', psutil.virtual_memory())
        yield ps
        print('mem2', psutil.virtual_memory())


for ps in iter_plot_todo():
    if ps.algo == 'global':
        continue
        plot_global(ps)
    else:
        plot_region(ps)
#    sys.exit(0)
