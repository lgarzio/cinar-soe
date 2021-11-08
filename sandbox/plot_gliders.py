#!/usr/bin/env python

"""
Author: Lori Garzio on 11/2/2021
Last modified: 11/2/2021
Glider data for OA State of the Ecosystem analysis
"""

import datetime as dt
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
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
import statistics
pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console
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


def main(sbuf, umf, ruf, lon_bounds, lat_bounds):
    bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'

    # for each glider, grab the data in the bottom 1m of the profile except when the profile is deeper than a
    # defined threshold (when the glider is off the shelf and no longer sampling the bottom)
    # we're assuming each profile samples the bottom, except off the shelf

    data = dict(pressure=np.array([]),
                lat=np.array([]),
                lon=np.array([]),
                pH=np.array([]))
    glider_data = dict(sbu=data,
                       ru=data)

    for rf in ruf:
        # RU30 - glider rated for 200m, so ignore profiles > 170m (presumably off the shelf)
        ru_data = glider_data['ru']
        ds = xr.open_dataset(rf)
        try:
            pid = np.unique(ds.profile_id.values)
            profile_identifier = 'profile_id'
            phvar = 'pH'
        except AttributeError:
            pid = np.unique(ds.profile_time.values)
            profile_identifier = 'profile_time'
            phvar = 'ph_total_shifted'
        for pd in pid:
            pidx = np.where(ds[profile_identifier].values == pd)[0]
            maxpress = np.nanmax(ds.pressure.values[pidx])
            if np.logical_and(maxpress > 30, maxpress < 170):
                idx = np.where(ds.pressure.values[pidx] > maxpress - 1)[0]
                ru_data['pressure'] = np.append(ru_data['pressure'], maxpress)
                ru_data['lat'] = np.append(ru_data['lat'], np.nanmean(ds.latitude.values[pidx][idx]))
                ru_data['lon'] = np.append(ru_data['lon'], np.nanmean(ds.longitude.values[pidx][idx]))
                ru_data['pH'] = np.append(ru_data['pH'], statistics.median(ds[phvar].values[pidx][idx]))

    for sf in sbuf:
        # SBU - glider rated for 350m but Charlie keeps it at 180m off the shelf, so ignore profiles > 170m
        sbu_data = glider_data['sbu']
        ds = xr.open_dataset(sf)
        ptime = np.unique(ds.profiletime)  # these data are actually grouped by yo (one down-up pair)
        for pt in ptime:
            dst = ds.sel(profiletime=pt)
            maxpress = np.nanmax(dst.pressure.values)
            if np.logical_and(maxpress > 30, maxpress < 170):
                idx = np.where(dst.pressure.values > maxpress-1)[0]
                sbu_data['pressure'] = np.append(sbu_data['pressure'], maxpress)
                sbu_data['lat'] = np.append(sbu_data['lat'], np.unique(dst.profilelat.values)[0])
                sbu_data['lon'] = np.append(sbu_data['lon'], np.unique(dst.profilelon.values)[0])
                sbu_data['pH'] = np.append(sbu_data['pH'], statistics.median(dst.ph_total.values[idx]))

    fig, ax = plt.subplots(figsize=(11, 8), subplot_kw=dict(projection=ccrs.Mercator()))
    plt.subplots_adjust(right=0.82)
    extent = [-78, -65, 35, 45]

    # add bathymetry
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))

    #levels = np.arange(-5000, 5100, 50)
    bath_lat = bathy.variables['lat'][:]
    bath_lon = bathy.variables['lon'][:]
    bath_elev = bathy.variables['elevation'][:]
    #plt.contourf(bath_lon, bath_lat, bath_elev,  levels, cmap=cmo.cm.topo, transform=ccrs.PlateCarree())

    levels = [-3000, -1000, -100]
    CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                     transform=ccrs.PlateCarree())
    ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

    add_map_features(ax, extent)

    sct = ax.scatter(glider_data['sbu']['lon'], glider_data['sbu']['lat'], c=glider_data['sbu']['pH'], marker='.',
                     s=75, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(), zorder=10)
    sct = ax.scatter(glider_data['ru']['lon'], glider_data['ru']['lat'], c=glider_data['ru']['pH'], marker='.',
                     s=75, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(), zorder=10)

    # Set colorbar height equal to plot height
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
    fig.add_axes(cax)

    # generate colorbar
    cb = plt.colorbar(sct, cax=cax)
    cb.set_label(label='Bottom pH (summer)')

    sfile = os.path.join('/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/bottom_ph_map_gliders.png')
    plt.savefig(sfile, dpi=200)
    plt.close()


if __name__ == '__main__':
    sbu_files = ['/Users/garzio/Documents/rucool/Saba/gliderdata/charlie/sbu01_04_final_May2021_calculated.nc',
                 '/Users/garzio/Documents/rucool/Saba/gliderdata/charlie/sbu01_05_final_July2021_calculated.nc']
    umaine_files = []
    ru_files = ['/Users/garzio/Documents/rucool/Saba/gliderdata/2019/ru30-20190717T1812/ru30-20190717T1812-delayed-dac.nc',
                '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210716T1804/delayed/ru30-20210716T1804-profile-sci-delayed_shifted.nc']
    lons = [-78, -65, -65, -78]
    lats = [35, 35, 45, 45]
    main(sbu_files, umaine_files, ru_files, lons, lats)
