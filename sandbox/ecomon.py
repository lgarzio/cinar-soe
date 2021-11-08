#!/usr/bin/env python

"""
Author: Lori Garzio on 10/18/2021
Last modified: 11/2/2021
CODAP-NA and glider datasets for OA analysis.
CODAP-NA dataset documented here: https://essd.copernicus.org/articles/13/2777/2021/
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
import PyCO2SYS as pyco2
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


def main(f, lon_bounds, lat_bounds, sbuf, umf, ruf, ecomon_files):
    bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'

    # get EcoMon data
    ecomon_data = dict(pressure=np.array([]),
                       lat=np.array([]),
                       lon=np.array([]),
                       pH=np.array([]))
    for ef in ecomon_files:
        df_ecomon = pd.read_csv(ef)
        df_ecomon.replace(-999, np.nan, inplace=True)
        df_ecomon.dropna(subset=['pH_TS_20C'], inplace=True)
        df_ecomon.dropna(subset=['Depth_Bottom_meters'], inplace=True)

        # combine Station_ID and Cast_number to get a unique profile
        df_ecomon['profile_id'] = df_ecomon.Station_ID.astype(str) + '_' + df_ecomon.Cast_number.astype(str)

        # for each profile, find the deepest sample and compare the sample depth to the water column depth
        # keep the sample if it's within 20% of the bottom depth
        profile_num = np.unique(df_ecomon.profile_id)
        for i, pn in enumerate(profile_num):
            df_pn = df_ecomon[df_ecomon.profile_id == pn]
            maxdepth = np.nanmax(df_pn.Depth_meters)
            if maxdepth > 10:  # sampling depth has to be >10m
                pn_coords = [df_pn[df_pn.Depth_meters == maxdepth].Longitude_Dec_Deg.values[0],
                             df_pn[df_pn.Depth_meters == maxdepth].Latitude_Dec_Deg.values[0]]
                bottom_depth = np.nanmax(df_pn.Depth_Bottom_meters)
                bottom_threshold = bottom_depth * .8
                if maxdepth > bottom_threshold:
                    ecomon_data['pressure'] = np.append(ecomon_data['pressure'], maxdepth)
                    ecomon_data['lat'] = np.append(ecomon_data['lat'], pn_coords[1])
                    ecomon_data['lon'] = np.append(ecomon_data['lon'], pn_coords[0])

                    # calculate pH at in situ temperature, pressure, and salinity
                    df_pn_row = df_pn[df_pn.Depth_meters == maxdepth]
                    ph_20c = df_pn_row.pH_TS_20C.values[0]
                    par1_type = 3
                    ta = df_pn_row['TA_umol/kg'].values[0]
                    par2_type = 1

                    kwargs = dict(salinity=df_pn_row.CTDSAL_PSS78.values[0],
                                  temperature=20,
                                  temperature_out=df_pn_row.CTDTEMP_ITS90.values[0],
                                  pressure=0,
                                  pressure_out=df_pn_row.CTDPRES_dbar.values[0],
                                  opt_pH_scale=1,
                                  opt_k_carbonic=4,
                                  opt_k_bisulfate=1,
                                  opt_total_borate=1,
                                  opt_k_fluoride=2)

                    results = pyco2.sys(ph_20c, ta, par1_type, par2_type, **kwargs)
                    ecomon_data['pH'] = np.append(ecomon_data['pH'], results['pH_out'])

    # plot
    fig, ax = plt.subplots(figsize=(11, 8), subplot_kw=dict(projection=ccrs.Mercator()))
    plt.subplots_adjust(right=0.82)

    # define bathymetry levels and data
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

    # plot ECOA data
    sct = ax.scatter(ecoa_data['lon'], ecoa_data['lat'], c=ecoa_data['pH'], marker='.',
                     s=100, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(), zorder=10)

    # Set colorbar height equal to plot height
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
    fig.add_axes(cax)

    # generate colorbar
    cb = plt.colorbar(sct, cax=cax)
    cb.set_label(label='Bottom pH (summer)')

    sfile = os.path.join(os.path.dirname(f), 'bottom_ph_map.png')
    plt.savefig(sfile, dpi=200)
    plt.close()


if __name__ == '__main__':
    fname = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/CODAP_NA_v2021.nc'
    lons = [-78, -65, -65, -78]
    lats = [35, 35, 45, 45]
    sbu_files = ['/Users/garzio/Documents/rucool/Saba/gliderdata/charlie/sbu01_04_final_May2021_calculated.nc',
                 '/Users/garzio/Documents/rucool/Saba/gliderdata/charlie/sbu01_05_final_July2021_calculated.nc']
    umaine_files = []
    ru_files = [
        '/Users/garzio/Documents/rucool/Saba/gliderdata/2019/ru30-20190717T1812/ru30-20190717T1812-delayed-dac.nc',
        '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210716T1804/delayed/ru30-20210716T1804-profile-sci-delayed_shifted.nc']
    ecomon = ['/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/EcoMon/33GG20190815-GU1902_data.csv',
              '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/EcoMon/33HH20190522-HB1902_data.csv']
    main(fname, lons, lats, sbu_files, umaine_files, ru_files, ecomon)
