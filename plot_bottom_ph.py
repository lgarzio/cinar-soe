#!/usr/bin/env python

"""
Author: Lori Garzio on 10/18/2021
Last modified: 11/16/2022
Plot summer bottom-water pH using CODAP-NA, EcoMon, and glider datasets.
CODAP-NA dataset documented here: https://essd.copernicus.org/articles/13/2777/2021/
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cmocean as cmo
import pickle
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console
plt.rcParams.update({'font.size': 13})


def add_map_features(axis, extent, edgecolor=None):
    edgecolor = edgecolor or 'black'

    land = cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor=edgecolor, facecolor='tan')

    state_lines = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none')

    # Axes properties and features
    axis.set_extent(extent)
    axis.add_feature(land)
    axis.add_feature(cfeature.RIVERS)
    axis.add_feature(cfeature.LAKES)
    axis.add_feature(cfeature.BORDERS)
    axis.add_feature(state_lines, zorder=11, edgecolor=edgecolor)

    # Gridlines and grid labels
    gl = axis.gridlines(
        draw_labels=True,
        linewidth=.5,
        color='black',
        alpha=0.25,
        linestyle='--'
    )

    gl.top_labels = gl.right_labels = False
    gl.xlabel_style = {'size': 11, 'color': 'black'}
    gl.ylabel_style = {'size': 11, 'color': 'black'}
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER


def main(cruise_file, glider_file, savedir):
    bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    extent = [-78, -65, 35, 45]
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))

    with open(cruise_file, 'rb') as handle:
        cdata = pickle.load(handle)

    with open(glider_file, 'rb') as handle:
        gdata = pickle.load(handle)

    # combine datasets
    lon_full = np.append(cdata['lon'], gdata['lon'])
    lat_full = np.append(cdata['lat'], gdata['lat'])
    data_full = np.append(cdata['pH'], gdata['pH'])

    # plot map of summer bottom pH data
    fig, ax = plt.subplots(figsize=(11, 8), subplot_kw=dict(projection=ccrs.Mercator()))
    plt.subplots_adjust(right=0.82)

    # define bathymetry levels and data
    bath_lat = bathy.variables['lat'][:]
    bath_lon = bathy.variables['lon'][:]
    bath_elev = bathy.variables['elevation'][:]

    levels = [-3000, -1000, -100]
    CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                     transform=ccrs.PlateCarree())
    ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

    add_map_features(ax, extent)

    sct = ax.scatter(lon_full, lat_full, c=data_full, marker='.', vmin=7.7, vmax=8.1,
                     s=100, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(), zorder=10)

    # Set colorbar height equal to plot height
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
    fig.add_axes(cax)

    # generate colorbar
    cb = plt.colorbar(sct, cax=cax, extend='both')
    cb.set_label(label='Bottom pH (summer)')

    sfile = os.path.join(savedir, 'bottom_pH_map_updated.png')
    plt.savefig(sfile, dpi=200)
    plt.close()


if __name__ == '__main__':
    cruise_data = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/pickle/cruise_summer_bottom_data.pickle'
    glider_data = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/pickle/glider_summer_bottom_data.pickle'
    save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data'
    main(cruise_data, glider_data, save_directory)
