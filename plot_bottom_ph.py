#!/usr/bin/env python

"""
Author: Lori Garzio on 10/18/2021
Last modified: 12/4/2023
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
import cmocean as cmo
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console
plt.rcParams.update({'font.size': 15})


def main(cruise_file, glider_file, ab, savedir):
    bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    extent = [-78, -65, 35, 45]
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))

    cdata = xr.open_dataset(cruise_file)
    gdata = xr.open_dataset(glider_file)

    # locations of "hot spot" boxes (Gulf of Mexico and NJ)
    gom_box_lons = [-70.4, -69.25, -69.7, -70.85, -70.4]
    gom_box_lats = [43.36, 43.1, 42.15, 42.4, 43.36]
    nj_box_lons = [-73.86, -72.75, -73.08, -74.2, -73.86]
    nj_box_lats = [40.4, 40.14, 39.28, 39.55, 40.4]

    # combine datasets
    lon_full = np.append(cdata.lon.values, gdata.lon.values)
    lat_full = np.append(cdata.lat.values, gdata.lat.values)
    data_full = np.append(cdata.pH_bottom.values, gdata.pH_bottom.values)
    time_full = np.append(cdata.time.values, gdata.time.values)
    years = np.unique(pd.to_datetime(time_full).year)
    min_year = np.min(years)
    max_year = np.max(years)
    year_list = [y for y in range(min_year, max_year + 1)]

    # plot map of summer bottom pH data for entire dataset
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

    cf.add_map_features(ax, extent)

    sct = ax.scatter(lon_full, lat_full, c=data_full, marker='.', vmin=7.7, vmax=8.1,
                     s=100, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(), zorder=10)

    # Set colorbar height equal to plot height
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
    fig.add_axes(cax)

    # generate colorbar
    cb = plt.colorbar(sct, cax=cax, extend='both')
    cb.set_label(label='Bottom pH (summer)')

    if ab:
        ax.fill(gom_box_lons, gom_box_lats, color='none', edgecolor='cyan', linewidth=2.5,
                transform=ccrs.PlateCarree(), zorder=10)
        ax.fill(nj_box_lons, nj_box_lats, color='none', edgecolor='cyan', linewidth=2.5,
                transform=ccrs.PlateCarree(), zorder=10)

        sfile = os.path.join(savedir, f'bottom_pH_map_{min_year}-{max_year}-boxes.png')
    else:
        sfile = os.path.join(savedir, f'bottom_pH_map_{min_year}-{max_year}.png')

    plt.savefig(sfile, dpi=200)
    plt.close()

    # plot maps of summer bottom pH data for each year
    for y in year_list:
        # grab data for the year
        idx = np.where(pd.to_datetime(time_full).year == y)[0]

        fig, ax = plt.subplots(figsize=(11, 8), subplot_kw=dict(projection=ccrs.Mercator()))
        plt.subplots_adjust(right=0.82)

        # add bathymetry
        CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                         transform=ccrs.PlateCarree())
        ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

        cf.add_map_features(ax, extent)

        sct = ax.scatter(lon_full[idx], lat_full[idx], c=data_full[idx], marker='.', vmin=7.7, vmax=8.1,
                         s=100, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(), zorder=10)

        # Set colorbar height equal to plot height
        divider = make_axes_locatable(ax)
        cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
        fig.add_axes(cax)

        # generate colorbar
        cb = plt.colorbar(sct, cax=cax, extend='both')
        cb.set_label(label='Bottom pH')

        # add title
        plt.suptitle(f'Summer {y}', x=.46, y=.93)

        if ab:
            ax.fill(gom_box_lons, gom_box_lats, color='none', edgecolor='cyan', linewidth=2.5,
                    transform=ccrs.PlateCarree(), zorder=10)
            ax.fill(nj_box_lons, nj_box_lats, color='none', edgecolor='cyan', linewidth=2.5,
                    transform=ccrs.PlateCarree(), zorder=10)

        if len(idx) > 0:
            versions = ['b', 'c', 'd', 'e', 'f', 'g']
        else:
            versions = ['b', 'c', 'd']

        for v in versions:
            if ab:
                sfile = os.path.join(savedir, 'bottom_pH_years_boxes', f'bottom_pH_map_{y}{v}-boxes.png')
            else:
                sfile = os.path.join(savedir, 'bottom_pH_years', f'bottom_pH_map_{y}{v}.png')
            plt.savefig(sfile, dpi=200)
        plt.close()


if __name__ == '__main__':
    cruise_data = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/output_nc/vessel_based_summer_bottom_OA_data.nc'
    glider_data = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/output_nc/glider_based_summer_bottom_OA_data.nc'
    add_boxes = False
    save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/plots2023'
    main(cruise_data, glider_data, add_boxes, save_directory)
