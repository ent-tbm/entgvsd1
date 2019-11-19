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
import psutil
import fcntl

"""OBSOLETE Python plotting"""

def plot_file(ifname, ofname):
    """region.center_ll:
        (lat, lon) of center of region to plot"""
    print('Plotting {}'.format(ifname))

    with netCDF4.Dataset(ifname) as nc:

        nlon = nc.variables['lon'].shape[0]
        nlat = nc.variables['lat'].shape[0]

        screenres = (1800,2880)    # y,x
        dotpitch = 1./220    # 220 dpi monitor

        # Read lons and lats (cell centers)
        lon = nc.variables['lon'][:]
        lat = nc.variables['lat'][:]

        # Get data variable name
        lonlat = {'lon', 'lat'}
        for vname,ncvar in nc.variables.items():
            if vname not in lonlat:
                break

        val = ncvar[:]
        val[val==0] = np.nan
        plotter = giss.plot.LonLatPlotter(lon, lat, boundaries=False)
     
        basemap = giss.basemap.global_map()

        # One plot per page (figsize must be in inches)
        figure = matplotlib.pyplot.figure(figsize=(screenres[1]*dotpitch, screenres[0]*dotpitch))

        ax = figure.add_subplot(111)

        cb_args=dict()
        if ifname.endswith('_lc_lr.nc'):
            plot_args=dict(vmin=0.,vmax=1.0)

        elif ifname.endswith('_lai_lr.nc'):
            plot_args=dict(vmin=0.,vmax=7.0)
            plot_args['cmap'] = giss.plot.read_cpt('NYT_drought14b.cpt').cmap
        else:
            raise ValueError('Unknown plot type (lc vs lai)')

        _title = nc.title if hasattr(nc,'title') else os.path.split(ifname)[1]
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
        print('    Writing {}'.format(ofname))
        figure.savefig('_tmp.png', dpi=1./dotpitch, transparent=False)
        os.rename('_tmp.png', ofname)   # Write atomically


def plot_all():
    for path,dirnames,filenames in os.walk('lc_lai_ent'):
        for fname in sorted(filenames):
            if fname.endswith('_lr.nc'):
                ifname = os.path.join(path,fname)
                ofname = '{}.png'.format(os.path.splitext(ifname)[0][:-3])
                print(ifname,ofname)

                if (not os.path.exists(ofname)) or (os.path.getmtime(ifname) > os.path.getmtime(ofname)):
                    # Check for ourselves if we are runnig low on memory.
                    # If it's this low, probably another plotregion is running.
                    # Exit with no error, so our caller doesn't keep calling us.
                    avail = psutil.virtual_memory().available
                    print('********** avail', avail)
                    if avail < 1000000000:
                        sys.exit(1)

                    print('mem1', psutil.virtual_memory())
                    plot_file(ifname, ofname)
                    print('mem2', psutil.virtual_memory())

lock_fname = 'plotregion.lock'
try:
    # http://tilde.town/~cristo/file-locking-in-python.html
    lock_file = open(lock_fname, 'w+')
    fcntl.flock(lock_file, fcntl.LOCK_EX | fcntl.LOCK_NB)
    plot_all()
except BlockingIOError:
    sys.stderr.write('Cannot open lockfile {}, quitting!\n'.format(lock_fname))
    sys.exit(0)    # Tell caller to quit
finally:
    fcntl.flock(lock_file, fcntl.LOCK_UN)
    os.remove(lock_fname)

