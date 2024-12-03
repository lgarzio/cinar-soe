#!/usr/bin/env python

"""
Author: Lori Garzio on 11/16/2022
Last modified: 11/12/2024
Plot in highlighted circles when summer bottom/surface aragonite saturation state (omega) is <= defined thresholds for
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
import cartopy.feature as cfeature
import cmocean as cmo
import functions.common as cf
import cool_maps.plot as cplt
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console
plt.rcParams.update({'font.size': 13})


def main(cruise_file, glider_file, ab, stype, savedir, thresh):
    savedir = os.path.join(savedir, stype)
    season_mapping = {'DJF': 'Winter',
                      'MAM': 'Spring',
                      'JJA': 'Summer',
                      'SON': 'Fall'}

    bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    extent = [-78, -65, 35, 45]
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))

    # define bathymetry levels and data
    bath_lat = bathy['lat'].values
    bath_lon = bathy['lon'].values
    bath_elev = bathy['elevation'].values
    levels = [-3000, -1000, -100]

    cdata = xr.open_dataset(cruise_file)
    gdata = xr.open_dataset(glider_file)

    # add month and season
    cdata['month'] = cdata['time.month']
    cdata['season'] = cdata['time.season']
    cdata['year'] = cdata['time.year']
    gdata['month'] = gdata['time.month']
    gdata['season'] = gdata['time.season']
    gdata['year'] = gdata['time.year']

    # convert to dataframe and merge
    cdf = cdata.to_dataframe()
    gdf = gdata.to_dataframe()
    cdf['deployment'] = ''
    gdf['cruise'] = ''
    cdf['sample_type'] = 'vessel'
    gdf['sample_type'] = 'glider'
    cdf = cdf[['cruise', 'deployment', 'lat', 'lon', 'depth', f'omega_{stype}', f'temperature_{stype}',
               'month', 'season', 'year', 'sample_type']]
    gdf = gdf[['cruise', 'deployment', 'lat', 'lon', 'depth', f'omega_{stype}', f'temperature_{stype}',
               'month', 'season', 'year', 'sample_type']]
    df = pd.concat([cdf, gdf])

    # get rid of nans
    df = df[~np.isnan(df[f'omega_{stype}'])]

    # # locations of "hot spot" boxes (Gulf of Maine and NJ)
    # gom_box_lons = [-70.4, -69.25, -69.7, -70.85, -70.4]
    # gom_box_lats = [43.36, 43.1, 42.15, 42.4, 43.36]
    # nj_box_lons = [-73.86, -72.75, -73.08, -74.2, -73.86]
    # nj_box_lats = [40.4, 40.14, 39.28, 39.55, 40.4]

    for season_code in ['JJA', 'MAM', 'SON', 'DJF']:
        print(season_code)
        season = season_mapping[season_code]
        df_season = df[df.season == season_code].copy()

        ###############################################################################################################
        # plot maps where values are below defined thresholds for select species
        for key, value in thresh.items():
            fig, ax = plt.subplots(figsize=(9, 8), subplot_kw=dict(projection=ccrs.Mercator()))
            plt.subplots_adjust(top=.91, bottom=0.08, right=.94, left=0.08)
            CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                             transform=ccrs.PlateCarree())
            ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

            #cplt.add_bathymetry(ax, bath_lon, bath_lat, bath_elev, method='topo_log', zorder=1)

            cf.add_map_features(ax, value['lims'], oceancolor=cfeature.COLORS['water'])

            # grab the data only within the defined limits for that species
            lon_bounds = [value['lims'][0], value['lims'][1], value['lims'][1], value['lims'][0]]
            lat_bounds = [value['lims'][2], value['lims'][2], value['lims'][3], value['lims'][3]]
            df_season['in_region'] = ''
            for i, row in df_season.iterrows():
                if Polygon(list(zip(lon_bounds, lat_bounds))).contains(Point(row.lon, row.lat)):
                    df_season.loc[i, 'in_region'] = 'yes'
                else:
                    df_season.loc[i, 'in_region'] = 'no'

            df_season_region = df_season[df_season['in_region'] == 'yes']

            # plot everything as empty circles
            sct = ax.scatter(df_season_region.lon, df_season_region.lat, c='None',
                             marker='o', edgecolor='lightgray', s=20, transform=ccrs.PlateCarree(), zorder=10)

            # plot the values less than the threshold and within the depth range as filled circles
            df_season_region_flag = df_season_region[df_season_region[f'omega_{stype}'] < value['omega_sensitivity']]
            df_season_region_flag = df_season_region_flag[(df_season_region_flag['depth'] >= value['depth_range'][0]) & (
                    df_season_region_flag['depth'] <= value['depth_range'][1])]

            sct = ax.scatter(df_season_region_flag.lon, df_season_region_flag.lat, c='darkcyan',
                             marker='o', s=20, transform=ccrs.PlateCarree(), zorder=10, label='<2023')

            # plot 2023 and 2024 in a different color
            colors = ['magenta', 'cyan']
            for i, yy in enumerate([2023, 2024]):
                # winter spans two years - grab the December data from the previous year
                if season_code == 'DJF':
                    # Dec from previous year
                    df_season_year1 = df_season_region_flag[(df_season_region_flag['year'] == (yy - 1)) & (df_season_region_flag['month'] == 12)]
                    # Jan and Feb for the current year
                    df_season_year2 = df_season_region_flag[
                        (df_season_region_flag['year'] == yy) & np.logical_or(df_season_region_flag['month'] == 1, df_season_region_flag['month'] == 2)]
                    df_season_year = pd.concat([df_season_year1, df_season_year2])
                else:
                    df_season_year = df_season_region_flag[(df_season_region_flag['year'] == yy) & (df_season_region_flag['season'] == season_code)]

                sct = ax.scatter(df_season_year.lon, df_season_year.lat, c=colors[i],
                                 marker='o', s=20, transform=ccrs.PlateCarree(), zorder=10, label=str(yy))

            plt.legend(loc='upper left')

            plt.title('{}: Omega below calcification sensitivity of {}\n{} (depth range: {}-{}m)'.format(
                season,
                value['omega_sensitivity'],
                value['long_name'],
                value['depth_range'][0],
                value['depth_range'][1]))

            min_year1 = np.nanmin(np.unique(df_season_region['year']))
            max_year1 = np.nanmax(np.unique(df_season_region['year']))
            sfile = os.path.join(savedir, f'{stype}_omega_map-{key}-{min_year1}-{max_year1}-{season.lower()}.png')
            plt.savefig(sfile, dpi=200)
            plt.close()

            # export a dataframe of the times the threshold was exceeded
            summary = pd.DataFrame(df_season_region_flag.index)
            df_days = summary.groupby(summary['time'].map(lambda x: x.day)).min()
            df_days.rename(columns={'time': 'date_threshold_reached'}, inplace=True)
            df_days.reset_index(inplace=True)
            df_days['date_threshold_reached'] = df_days['date_threshold_reached'].map(lambda t: t.strftime('%Y-%m-%d'))
            df_days.sort_values(by='date_threshold_reached', inplace=True)
            df_days.drop(columns=['time'], inplace=True)

            df_days.to_csv(os.path.join(savedir, f'{stype}_omega_days_threshold_reached-{key}-{season.lower()}.csv'), index=False)

            # make a plot for each year
            region_years = np.unique(df_season_region.year)
            for yy in region_years:
                fig, ax = plt.subplots(figsize=(9, 8), subplot_kw=dict(projection=ccrs.Mercator()))
                plt.subplots_adjust(top=.9, bottom=0.08, right=.94, left=0.08)
                CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                                 transform=ccrs.PlateCarree())
                ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

                #cplt.add_bathymetry(ax, bath_lon, bath_lat, bath_elev, method='topo_log', zorder=1)

                cf.add_map_features(ax, value['lims'], oceancolor=cfeature.COLORS['water'])

                # grab data for the year
                # winter spans two years - grab the December data from the previous year
                if season_code == 'DJF':
                    # Dec from previous year
                    df_season_year1 = df_season_region[(df_season_region['year'] == (yy - 1)) & (df_season_region['month'] == 12)]
                    # Jan and Feb for the current year
                    df_season_year2 = df_season_region[
                        (df_season_region['year'] == yy) & np.logical_or(df_season_region['month'] == 1, df_season_region['month'] == 2)]
                    df_season_year = pd.concat([df_season_year1, df_season_year2])
                else:
                    df_season_year = df_season_region[
                        (df_season_region['year'] == yy) & (df_season_region['season'] == season_code)]

                # plot everything as empty circles
                if len(df_season_year) > 0:
                    sct = ax.scatter(df_season_year.lon, df_season_year.lat,
                                     marker='o', c='lightgray', s=20, transform=ccrs.PlateCarree(), zorder=10)

                    c = 'darkcyan'
                    if yy == 2023:
                        c = 'magenta'
                    if yy == 2024:
                        c = 'cyan'

                    # plot the values less than the threshold and within the depth range as filled circles
                    df_season_year_flag = df_season_year[
                        df_season_year[f'omega_{stype}'] < value['omega_sensitivity']]
                    df_season_year_flag = df_season_year_flag[
                        (df_season_year_flag['depth'] >= value['depth_range'][0]) & (df_season_year_flag['depth'] <= value['depth_range'][1])]

                    sct = ax.scatter(df_season_year_flag.lon, df_season_year_flag.lat, c=c,
                                     marker='o', s=20, transform=ccrs.PlateCarree(), zorder=10)

                    plt.title('{} {}: Omega below calcification sensitivity of {}:\n{} (depth range: {}-{}m)'.format(
                        season,
                        yy,
                        value['omega_sensitivity'],
                        value['long_name'],
                        value['depth_range'][0],
                        value['depth_range'][1]))

                    sdir = os.path.join(savedir, 'years', key)
                    os.makedirs(sdir, exist_ok=True)
                    sfile = os.path.join(sdir, f'{stype}_omega_map-{key}-{season.lower()}{yy}.png')
                    plt.savefig(sfile, dpi=200)
                    plt.close()


if __name__ == '__main__':
    #cruise_data = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/output_nc/vessel_based_bottom_OA_data_2007_2023.nc'  # bottom
    #glider_data = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/output_nc/glider_based_bottom_OA_data_2019_2024.nc'  # bottom
    cruise_data = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/output_nc/vessel_based_surface_OA_data_2004_2023.nc'  # surface
    glider_data = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/output_nc/glider_based_surface_OA_data_2019_2024.nc'  # surface
    add_boxes = False
    stype = 'surface'  # bottom surface
    save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/plots2025/species_omega_sensitivity'
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

    main(cruise_data, glider_data, add_boxes, stype, save_directory, thresholds)
