import netCDF4
import giss.basemap
#import giss.modele
import giss.plot
import matplotlib.pyplot
import mpl_toolkits.basemap
import sys
import numpy as np
import os
import collections

Region = collections.namedtuple('Region',
    ('name',        # Descriptive name of the region
    'shortname',    # Used in filenames
    'center_ll'))   # (lat, lon) of the center of the region


def plot_file(ifname, ofname, region):
    """region.center_ll:
        (lat, lon) of center of region to plot"""
    print('Plotting {} ({})'.format(ifname, region.name))

    with netCDF4.Dataset(ifname) as nc:

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
        lonlat = {'lon', 'lat'}
        for vname,ncvar in nc.variables.items():
            if vname not in lonlat:
                break

        val = ncvar[_base[0]:_top[0],  _base[1]:_top[1]]
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

        # Plot multiple plots on one page
        # (figsize must be in inches)
        figure = matplotlib.pyplot.figure(figsize=(screenres[1]*dotpitch, screenres[0]*dotpitch))

        ax = figure.add_subplot(111)


        cb_args=dict()
        if ifname.endswith('_lc.nc'):
            plot_args=dict(vmin=0.,vmax=1.0)
        elif ifname.endswith('_lai.nc'):
            plot_args=dict(vmin=0.,vmax=7.0)
        else:
            raise ValueError('Unknown plot type (lc vs lai)')

        _title = nc.title if hasattr(nc,'title') else os.path.split(ifname)[1]
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
        print('    Writing {}'.format(ofname))
        figure.savefig(ofname, dpi=1./dotpitch, transparent=False)


# Convert (deg,min,sec) into decimal lat/lon
def _N(deg,min,sec):
    return deg + min/60. + sec/3600.
def _S(deg,min,sec):
    return -(deg + min/60. + sec/3600.)
def _E(deg,min,sec):
    return deg + min/60. + sec/3600.
def _W(deg,min,sec):
    return - (deg + min/60. + sec/3600.)

regions = (
    Region('Kansas', 'ks', (39., -98.)),    # (lat,lon)
    Region('Sub-Saharan Africa', 'afk', (_N(11,16,00.9), _E(18,59,44.2))),
    Region('Europe crop belt', 'eur', (_N(50,35,58.9), _E(19,50,00.1))),
    Region('Brazil crops/shrubland', 'bra', (_S(17,25,34.2), _W(53,0,40.5))),
)

for path,dirnames,filenames in os.walk('lc_lai_ent'):
    for fname in sorted(filenames):
        if fname.endswith('.nc'):
            for region in regions:
                ifname = os.path.join(path,fname)
                ofname = '{}_{}.png'.format(os.path.splitext(ifname)[0], region.shortname)
                if (not os.path.exists(ofname)) or (os.path.getmtime(ifname) > os.path.getmtime(ofname)):
                    plot_file(ifname, ofname, region)

#for region in regions:
#
#plot_file('lc_lai_ent/EntMM_lc_laimax_1kmx1km/11_c3_grass_per_lc.nc')
#plot_file('lc_lai_ent/EntMM_lc_laimax_1kmx1km/16_crops_c4_herb_lc.nc')
