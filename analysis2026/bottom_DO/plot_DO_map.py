#!/usr/bin/env python

"""
Author: Lori Garzio on 10/29/2025
Last modified: 10/29/2025
Quick plot of bottom DO maps
"""

import os
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cool_maps.plot as cplt
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console
plt.rcParams.update({'font.size': 15})


def main(f):
    sfile = os.path.join(os.path.dirname(f), 'glider_bottom_DO_map_2025.png')
    extent = [-76, -69, 37, 42.5]

    ds = xr.open_dataset(f)

    kwargs = dict()
    kwargs['figsize'] = (9, 8)
    kwargs['coast'] = 'full'  # low full
    kwargs['oceancolor'] = 'none'
    kwargs['decimal_degrees'] = True
    kwargs['bathymetry'] = True
    kwargs['bathymetry_method'] = 'topo_log'
    fig, ax = cplt.create(extent, **kwargs)

    ax.scatter(ds['lon'], ds['lat'], color='magenta', marker='.', s=6, transform=ccrs.PlateCarree(), zorder=30)

    plt.savefig(sfile, dpi=200)
    plt.close()


if __name__ == '__main__':
    file = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/dissolved_oxygen_gliders/glider_based_bottom_DO_data_2025.nc'
    main(file)
