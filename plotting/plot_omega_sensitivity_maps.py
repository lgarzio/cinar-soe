#!/usr/bin/env python

"""
Author: Lori Garzio on 11/16/2022
Last modified: 10/2/2026
Plot in highlighted circles when summer bottom/surface aragonite saturation state (omega) is <= defined thresholds for
key Mid-Atlantic species using CODAP-NA, EcoMon, and glider datasets.
CODAP-NA dataset documented here: https://essd.copernicus.org/articles/13/2777/2021/
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
import yaml
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


def main(species_list, stype, savedir):
    home_dir = os.getenv('HOME')
    savedir = os.path.join(home_dir, savedir, stype)
    os.makedirs(savedir, exist_ok=True)
    season_mapping = {'DJF': 'Winter',
                      'MAM': 'Spring',
                      'JJA': 'Summer',
                      'SON': 'Fall'}

    bathymetry = os.path.join(home_dir, 'rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc')
    extent = [-78, -65, 35, 45]
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))

    # define bathymetry levels and data
    bath_lat = bathy['lat'].values
    bath_lon = bathy['lon'].values
    bath_elev = bathy['elevation'].values
    levels = [-3000, -1000, -100]

    # grab the dataset from ERDDAP
    ru_server = 'https://rucool-sampling.marine.rutgers.edu/erddap'
    dataset_id = 'vessel_glider_surface_bottom_OA_data_NES'
    ds = cf.return_erddap_nc(ru_server, dataset_id)
    ds = ds.swap_dims({'row': 'time'})

    # add year, month, season
    ds['year'] = ds['time.year']
    ds['month'] = ds['time.month']
    ds['season'] = ds['time.season']

    # convert to dataframe and get rid of nans
    df = ds.to_dataframe().reset_index()
    df = df[['time', 'year', 'month', 'season', 'latitude', 'longitude', 'depth', f'temperature_{stype}', f'omega_{stype}']]
    df = df[~np.isnan(df[f'omega_{stype}'])]

    # get threshold configuration file
    cfile = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'configs', 'species_thresholds.yml')
    with open(cfile, "r") as file:
        thresh = yaml.safe_load(file)

    # arguments for the map for all plots
    kwargs = dict()
    kwargs['figsize'] = (9, 8)
    kwargs['coast'] = 'high'  # low high full
    kwargs['oceancolor'] = cfeature.COLORS['water']
    kwargs['decimal_degrees'] = True
    kwargs['bathymetry'] = False
    kwargs['padding'] = 0

    for season_code in ['JJA', 'MAM', 'SON', 'DJF']:
        print(season_code)
        season = season_mapping[season_code]
        df_season = df[df.season == season_code].copy()

        ###############################################################################################################
        # plot maps where values are below defined thresholds for select species
        for sp in species_list:
            values = thresh[sp]
            if stype not in values['plttype']:
                continue
            if values['lims'][0] == -78:
                region = 'mab'
            elif values['lims'][0] == -72:
                region = 'gom'
            fig, ax = cplt.create(extent, **kwargs)
            plt.subplots_adjust(top=.91, bottom=0.08, right=.94, left=0.08)
            CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                             transform=ccrs.PlateCarree())
            ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

            # grab the data only within the defined limits for that species
            lon_bounds = [values['lims'][0], values['lims'][1], values['lims'][1], values['lims'][0]]
            lat_bounds = [values['lims'][2], values['lims'][2], values['lims'][3], values['lims'][3]]
            df_season['in_region'] = ''
            for i, row in df_season.iterrows():
                if Polygon(list(zip(lon_bounds, lat_bounds))).contains(Point(row.longitude, row.latitude)):
                    df_season.loc[i, 'in_region'] = 'yes'
                else:
                    df_season.loc[i, 'in_region'] = 'no'

            df_season_region = df_season[df_season['in_region'] == 'yes']

            # remove points outside of the species depth range
            df_season_region_depth = df_season_region[(df_season_region['depth'] >= values['depth_range'][0]) & (
                    df_season_region['depth'] <= values['depth_range'][1])]

            # plot everything within the species depth range as empty circles
            sct = ax.scatter(df_season_region_depth.longitude, df_season_region_depth.latitude, c='None',
                             marker='o', edgecolor='lightgray', s=20, transform=ccrs.PlateCarree(), zorder=10)

            # plot the values less than the threshold as filled circles
            df_season_region_flag = df_season_region_depth[df_season_region_depth[f'omega_{stype}'] < values['omega_sensitivity']]
            # df_season_region_flag = df_season_region_flag[(df_season_region_flag['depth'] >= values['depth_range'][0]) & (
            #         df_season_region_flag['depth'] <= values['depth_range'][1])]

            sct = ax.scatter(df_season_region_flag.longitude, df_season_region_flag.latitude, c='darkcyan',
                             marker='o', s=20, transform=ccrs.PlateCarree(), zorder=10, label='<2023')

            # plot 2023, 2024, 2025 in a different color
            colors = ['magenta', 'cyan', '#ff7240']  # lime green #89F336, neon orange #ff7240
            highlightyrs = [2023, 2024, 2025]
            loopyrs = np.intersect1d(np.unique(df_season_region_depth.year), highlightyrs)

            for i, yy in enumerate(loopyrs):
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

                if yy == 2025:
                    if season_code != 'SON':
                        sct = ax.scatter(df_season_year.longitude, df_season_year.latitude, c=colors[i],
                                         marker='o', s=20, transform=ccrs.PlateCarree(), zorder=10, label=str(yy))
                else:
                    sct = ax.scatter(df_season_year.longitude, df_season_year.latitude, c=colors[i], marker='o', s=20, transform=ccrs.PlateCarree(), zorder=10, label=str(yy))

            plt.legend(loc='upper left', framealpha=1).set_zorder(20)

            plt.title('{}: Omega below calcification sensitivity of {}\n{} (depth range: {}-{}m)'.format(
                season,
                values['omega_sensitivity'],
                values['long_name'],
                values['depth_range'][0],
                values['depth_range'][1]))

            min_year1 = np.nanmin(np.unique(df_season_region_depth['year']))
            max_year1 = np.nanmax(np.unique(df_season_region_depth['year']))
            sfile = os.path.join(savedir, f'{stype}_omega_map-{sp}-{min_year1}-{max_year1}-{region}-{season.lower()}.png')
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

            df_days.to_csv(os.path.join(savedir, f'{stype}_omega_days_threshold_reached-{sp}-{region}-{season.lower()}.csv'), index=False)

            # make a plot for each year
            region_years = np.unique(df_season_region_depth.year)
            for yy in region_years:
                fig, ax = cplt.create(extent, **kwargs)
                plt.subplots_adjust(top=.9, bottom=0.08, right=.94, left=0.08)
                CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                                 transform=ccrs.PlateCarree())
                ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')
            
                # grab data for the year
                # winter spans two years - grab the December data from the previous year
                if season_code == 'DJF':
                    # Dec from previous year
                    df_season_year1 = df_season_region_depth[(df_season_region_depth['year'] == (yy - 1)) & (df_season_region_depth['month'] == 12)]
                    # Jan and Feb for the current year
                    df_season_year2 = df_season_region_depth[
                        (df_season_region_depth['year'] == yy) & np.logical_or(df_season_region_depth['month'] == 1, df_season_region_depth['month'] == 2)]
                    df_season_year = pd.concat([df_season_year1, df_season_year2])
                else:
                    df_season_year = df_season_region_depth[
                        (df_season_region_depth['year'] == yy) & (df_season_region_depth['season'] == season_code)]
            
                # plot everything as empty circles
                if len(df_season_year) > 0:
                    sct = ax.scatter(df_season_year.longitude, df_season_year.latitude,
                                     marker='o', c='lightgray', s=20, transform=ccrs.PlateCarree(), zorder=10)
            
                    c = 'darkcyan'
                    if yy == 2023:
                        c = 'magenta'
                    if yy == 2024:
                        c = 'cyan'
                    if yy == 2025:
                        c = '#ff7240'
            
                    # plot the values less than the threshold and within the depth range as filled circles
                    df_season_year_flag = df_season_year[
                        df_season_year[f'omega_{stype}'] < values['omega_sensitivity']]
                    df_season_year_flag = df_season_year_flag[
                        (df_season_year_flag['depth'] >= values['depth_range'][0]) & (df_season_year_flag['depth'] <= values['depth_range'][1])]
            
                    sct = ax.scatter(df_season_year_flag.longitude, df_season_year_flag.latitude, c=c,
                                     marker='o', s=20, transform=ccrs.PlateCarree(), zorder=10)
            
                    plt.title('{} {}: Omega below calcification sensitivity of {}:\n{} (depth range: {}-{}m)'.format(
                        season,
                        yy,
                        values['omega_sensitivity'],
                        values['long_name'],
                        values['depth_range'][0],
                        values['depth_range'][1]))
            
                    sdir = os.path.join(savedir, 'years', sp)
                    os.makedirs(sdir, exist_ok=True)
                    sfile = os.path.join(sdir, f'{stype}_omega_map-{sp}-{season.lower()}{yy}.png')
                    plt.savefig(sfile, dpi=200)
                    plt.close()


if __name__ == '__main__':
    species = ['atlantic_sea_scallop'] # ['atlantic_sea_scallop', 'longfin_squid', 'lobster', 'cod', 'pteropod']
    stype = 'bottom'  # bottom surface
    save_directory = 'rucool/Saba/NOAA_SOE/data/plots2027/species_omega_sensitivity'
    main(species, stype, save_directory)
