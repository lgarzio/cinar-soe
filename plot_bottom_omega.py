#!/usr/bin/env python

"""
Author: Lori Garzio on 11/16/2022
Last modified: 11/16/2022
Plot summer bottom-water aragonite saturation state (omega) using CODAP-NA, EcoMon, and glider datasets.
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
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
import statistics
import copy
import PyCO2SYS as pyco2
import functions.common as cf
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


def main(cruise_file, glider_file, savedir, thresh):
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
    data_full = np.append(cdata['aragonite'], gdata['aragonite'])
    depth_full = np.append(cdata['depth'], gdata['depth'])

    # get rid of nans
    nonan_idx = np.where(~np.isnan(data_full))[0]
    lon_full = lon_full[nonan_idx]
    lat_full = lat_full[nonan_idx]
    data_full = data_full[nonan_idx]
    depth_full = depth_full[nonan_idx]

    # plot map of summer bottom omega data
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

    sct = ax.scatter(lon_full, lat_full, c=data_full, marker='.', vmin=.75, vmax=2.25,
                     s=50, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(), zorder=10)

    # Set colorbar height equal to plot height
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
    fig.add_axes(cax)

    # generate colorbar
    cb = plt.colorbar(sct, cax=cax, extend='both')
    cb.set_label(label='Bottom Aragonite Saturation State (summer)')

    sfile = os.path.join(savedir, 'bottom_omega_map.png')
    plt.savefig(sfile, dpi=200)
    plt.close()

    print('Cruise data min: {}'.format(np.nanmin(cdata['aragonite'])))
    print('Cruise data max: {}'.format(np.nanmax(cdata['aragonite'])))
    print('Glider data min: {}'.format(np.nanmin(gdata['aragonite'])))
    print('Glider data max: {}'.format(np.nanmax(gdata['aragonite'])))

    # plot maps where values are below defined thresholds for select species
    for key, value in thresh.items():
        fig, ax = plt.subplots(figsize=(11, 8), subplot_kw=dict(projection=ccrs.Mercator()))
        plt.subplots_adjust(right=0.82)
        CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                         transform=ccrs.PlateCarree())
        ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

        add_map_features(ax, extent)

        # plot the values above the threshold and outside of the organism's depth range as empty circles
        thresh_idx = np.where(data_full > value['omega_sensitivity'])[0]
        depth_idx = np.where(np.logical_or(depth_full < value['depth_range'][0], depth_full > value['depth_range'][1]))[0]

        empty_idx = np.unique(np.append(thresh_idx, depth_idx))

        sct = ax.scatter(lon_full[empty_idx], lat_full[empty_idx], c='None',
                         marker='o', edgecolor='lightgray', s=20, transform=ccrs.PlateCarree(), zorder=10)

        # plot the values less than the threshold and within the depth range as filled circles
        mask = np.ones(len(lon_full), bool)
        mask[empty_idx] = 0

        sct = ax.scatter(lon_full[mask], lat_full[mask], c='darkcyan',
                         marker='o', s=20, transform=ccrs.PlateCarree(), zorder=10)

        plt.title('Omega below calcification sensitivity of {}:\n{} (depth range: {}-{}m)'.format(value['omega_sensitivity'],
                                                                                                  value['long_name'],
                                                                                                  value['depth_range'][0],
                                                                                                  value['depth_range'][1]))

        sfile = os.path.join(savedir, 'bottom_omega_map-{}.png'.format(key))
        plt.savefig(sfile, dpi=200)
        plt.close()


if __name__ == '__main__':
    cruise_data = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/pickle/cruise_summer_bottom_data.pickle'
    glider_data = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/pickle/glider_summer_bottom_data.pickle'
    save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data'
    thresholds = dict(
        atlantic_sea_scallop=dict(
            long_name='Atlantic Sea Scallop',
            omega_sensitivity=1.1,
            depth_range=[30, 90]
        ),
        lobster=dict(
            long_name='American Lobster',
            omega_sensitivity=1.45,
            depth_range=[4, 50]
        ),
        cod=dict(
            long_name='Atlantic Cod',
            omega_sensitivity=1,
            depth_range=[0, 500]
        )
    )

    main(cruise_data, glider_data, save_directory, thresholds)
