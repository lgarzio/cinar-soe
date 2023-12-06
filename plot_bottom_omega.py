#!/usr/bin/env python

"""
Author: Lori Garzio on 11/16/2022
Last modified: 12/5/2023
Plot in highlighted circles when summer bottom-water aragonite saturation state (omega) is <= defined thresholds for
key Mid-Atlantic species using CODAP-NA, EcoMon, and glider datasets.
CODAP-NA dataset documented here: https://essd.copernicus.org/articles/13/2777/2021/
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs
import cmocean as cmo
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console
plt.rcParams.update({'font.size': 13})


def main(cruise_file, glider_file, ab, savedir, thresh):
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
    data_full = np.append(cdata.omega_bottom.values, gdata.omega_bottom.values)
    depth_full = np.append(cdata.depth.values, gdata.depth.values)
    time_full = np.append(cdata.time.values, gdata.time.values)
    years = np.unique(pd.to_datetime(time_full).year)
    min_year = np.min(years)
    max_year = np.max(years)
    year_list = [y for y in range(min_year, max_year + 1)]

    # get rid of nans
    nonan_idx = np.where(~np.isnan(data_full))[0]
    lon_full = lon_full[nonan_idx]
    lat_full = lat_full[nonan_idx]
    data_full = data_full[nonan_idx]
    depth_full = depth_full[nonan_idx]
    time_full = time_full[nonan_idx]

    # plot map of summer bottom omega for the entire dataset
    fig, ax = plt.subplots(figsize=(9, 8), subplot_kw=dict(projection=ccrs.Mercator()))
    plt.subplots_adjust(top=.92, bottom=0.08, right=.96, left=0)

    # define bathymetry levels and data
    bath_lat = bathy.variables['lat'][:]
    bath_lon = bathy.variables['lon'][:]
    bath_elev = bathy.variables['elevation'][:]

    levels = [-3000, -1000, -100]
    CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                     transform=ccrs.PlateCarree())
    ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

    cf.add_map_features(ax, extent)

    sct = ax.scatter(lon_full, lat_full, c=data_full, marker='.', vmin=.75, vmax=2.25,
                     s=100, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(), zorder=10)

    # Set colorbar height equal to plot height
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
    fig.add_axes(cax)

    # generate colorbar
    cb = plt.colorbar(sct, cax=cax, extend='both')
    cb.set_label(label='Bottom Aragonite Saturation State (summer)')

    if ab:
        ax.fill(gom_box_lons, gom_box_lats, color='none', edgecolor='cyan', linewidth=2.5,
                transform=ccrs.PlateCarree(), zorder=10)
        ax.fill(nj_box_lons, nj_box_lats, color='none', edgecolor='cyan', linewidth=2.5,
                transform=ccrs.PlateCarree(), zorder=10)

        sfile = os.path.join(savedir, f'bottom_omega_map_{min_year}-{max_year}-boxes.png')
    else:
        sfile = os.path.join(savedir, f'bottom_omega_map_{min_year}-{max_year}.png')

    plt.savefig(sfile, dpi=200)
    plt.close()

    print('Cruise data min: {}'.format(np.nanmin(cdata['omega_bottom'])))
    print('Cruise data max: {}'.format(np.nanmax(cdata['omega_bottom'])))
    print('Glider data min: {}'.format(np.nanmin(gdata['omega_bottom'])))
    print('Glider data max: {}'.format(np.nanmax(gdata['omega_bottom'])))

    ###################################################################################################################

    # plot maps of summer bottom omega data for each year
    for y in year_list:
        # grab data for the year
        idx = np.where(pd.to_datetime(time_full).year == y)[0]

        fig, ax = plt.subplots(figsize=(9, 8), subplot_kw=dict(projection=ccrs.Mercator()))
        plt.subplots_adjust(top=.92, bottom=0.08, right=.96, left=0)

        # add bathymetry
        CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                         transform=ccrs.PlateCarree())
        ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

        cf.add_map_features(ax, extent)

        sct = ax.scatter(lon_full[idx], lat_full[idx], c=data_full[idx], marker='.', vmin=.75, vmax=2.25,
                         s=100, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(), zorder=10)

        # Set colorbar height equal to plot height
        divider = make_axes_locatable(ax)
        cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
        fig.add_axes(cax)

        # generate colorbar
        cb = plt.colorbar(sct, cax=cax, extend='both')
        cb.set_label(label='Bottom Aragonite Saturation State')

        # add title
        plt.suptitle(f'Summer {y}', x=.48, y=.96)

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
                sfile = os.path.join(savedir, 'bottom_omega_years_boxes', f'bottom_omega_map_{y}{v}-boxes.png')
            else:
                sfile = os.path.join(savedir, 'bottom_omega_years', f'bottom_omega_map_{y}{v}.png')
            plt.savefig(sfile, dpi=200)
        plt.close()

    ###################################################################################################################
    # plot maps where values are below defined thresholds for select species
    savedir = os.path.join(savedir, 'species_omega_thresholds')
    for key, value in thresh.items():
        fig, ax = plt.subplots(figsize=(9, 8), subplot_kw=dict(projection=ccrs.Mercator()))
        plt.subplots_adjust(top=.92, bottom=0.08, right=.94, left=0.08)
        CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                         transform=ccrs.PlateCarree())
        ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

        cf.add_map_features(ax, value['lims'])

        # grab the data only within the defined limits for that species
        lon_bounds = [value['lims'][0], value['lims'][1], value['lims'][1], value['lims'][0]]
        lat_bounds = [value['lims'][2], value['lims'][2], value['lims'][3], value['lims'][3]]
        region_idx = []
        for i, lon in enumerate(lon_full):
            if Polygon(list(zip(lon_bounds, lat_bounds))).contains(Point(lon, lat_full[i])):
                region_idx.append(i)

        lon_full_region = lon_full[region_idx]
        lat_full_region = lat_full[region_idx]
        data_full_region = data_full[region_idx]
        depth_full_region = depth_full[region_idx]
        time_full_region = time_full[region_idx]

        # plot the values above the threshold and outside of the organism's depth range as empty circles
        thresh_idx = np.where(data_full_region > value['omega_sensitivity'])[0]
        depth_idx = np.where(np.logical_or(depth_full_region < value['depth_range'][0], depth_full_region > value['depth_range'][1]))[0]

        empty_idx = np.unique(np.append(thresh_idx, depth_idx))

        sct = ax.scatter(lon_full_region[empty_idx], lat_full_region[empty_idx], c='None',
                         marker='o', edgecolor='lightgray', s=20, transform=ccrs.PlateCarree(), zorder=10)

        # plot the values less than the threshold and within the depth range as filled circles
        mask = np.ones(len(lon_full_region), bool)
        mask[empty_idx] = 0

        sct = ax.scatter(lon_full_region[mask], lat_full_region[mask], c='darkcyan',
                         marker='o', s=20, transform=ccrs.PlateCarree(), zorder=10)

        plt.title('Omega below calcification sensitivity of {}:\n{} (depth range: {}-{}m)'.format(value['omega_sensitivity'],
                                                                                                  value['long_name'],
                                                                                                  value['depth_range'][0],
                                                                                                  value['depth_range'][1]))

        sfile = os.path.join(savedir, f'bottom_omega_map-{key}-{min_year}-{max_year}.png')
        plt.savefig(sfile, dpi=200)
        plt.close()

        # export a dataframe of the times the threshold was exceeded
        df = pd.DataFrame(sorted(time_full_region[mask]), columns=['time'])
        df_days = df.groupby(df['time'].map(lambda x: x.day)).min()
        df_days.rename(columns={'time': 'date_threshold_reached'}, inplace=True)
        df_days.reset_index(inplace=True)
        df_days['date_threshold_reached'] = df_days['date_threshold_reached'].map(lambda t: t.strftime('%Y-%m-%d'))
        df_days.sort_values(by='date_threshold_reached', inplace=True)
        df_days.drop(columns=['time'], inplace=True)

        df['time'] = df['time'].map(lambda t: t.strftime('%Y-%m-%d %H:%M:%S'))

        df.to_csv(os.path.join(savedir, 'bottom_omega_full-{}.csv'.format(key)), index=False)
        df_days.to_csv(os.path.join(savedir, 'bottom_omega_days_threshold_reached-{}.csv'.format(key)), index=False)


if __name__ == '__main__':
    cruise_data = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/output_nc/vessel_based_summer_bottom_OA_data.nc'
    glider_data = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/output_nc/glider_based_summer_bottom_OA_data.nc'
    add_boxes = False
    save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/plots2024'
    thresholds = dict(
        atlantic_sea_scallop=dict(
            long_name='Atlantic Sea Scallop',
            omega_sensitivity=1.1,
            depth_range=[25, 200],
            lims=[-78, -68, 36, 42]
        ),
        longfin_squid=dict(
            long_name='Longfin Squid',
            omega_sensitivity=0.96,
            depth_range=[0, 400],
            lims=[-78, -68, 36, 42]
        ),
        lobster=dict(
            long_name='American Lobster',
            omega_sensitivity=1.09,
            depth_range=[10, 700],
            lims=[-72, -65, 40.5, 45]
        ),
        cod=dict(
            long_name='Atlantic Cod',
            omega_sensitivity=1.19,
            depth_range=[10, 200],
            lims=[-72, -65, 40.5, 45]
        )
    )

    main(cruise_data, glider_data, add_boxes, save_directory, thresholds)
