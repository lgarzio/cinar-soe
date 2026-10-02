#!/usr/bin/env python

"""
Author: Lori Garzio on 10/18/2021
Last modified: 10/2/2026
Plot bottom/surface omega/pH maps using CODAP-NA, EcoMon, and glider datasets.
ERDDAP server: https://rucool-sampling.marine.rutgers.edu/erddap/tabledap/vessel_glider_surface_bottom_OA_data_NES.html
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cmocean as cmo
import functions.common as cf
import cool_maps.plot as cplt
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console
plt.rcParams.update({'font.size': 15})


def main(stype, variable, vers, clims, addns, savedir):
    home_dir = os.getenv('HOME')
    savedir = os.path.join(home_dir, savedir, f'{stype}_{variable}')
    os.makedirs(savedir, exist_ok=True)
    season_mapping = {'DJF': 'Winter',
                      'MAM': 'Spring',
                      'JJA': 'Summer',
                      'SON': 'Fall'}

    # make a summary file
    rows = []
    columns = ['year', 'season', 'glider_date_range', 'vessel_date_range', f'min_glider_{variable}_{stype}', f'max_glider_{variable}_{stype}', f'min_glider_temp_{stype}',
               f'max_glider_temp_{stype}', 'n_glider', f'min_vessel_{variable}_{stype}', f'max_vessel_{variable}_{stype}',
               f'min_vessel_temp_{stype}', f'max_vessel_temp_{stype}', 'n_vessel', 'cruises', 'deployments']

    # for bathymetry contours
    bathymetry = f'{home_dir}/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
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

    # arguments for the map for all plots
    kwargs = dict()
    kwargs['figsize'] = (9, 8)
    kwargs['coast'] = 'high'  # low high full
    kwargs['oceancolor'] = cfeature.COLORS['water']
    kwargs['decimal_degrees'] = True
    kwargs['bathymetry'] = False
    kwargs['padding'] = 0
    
    for season_code in np.unique(ds.season):
        print(season_code)
        season = season_mapping[season_code]
        ds_season = ds.where(ds["season"] == season_code, drop=True)

        if vers == 'glider_only':
            ds_season = ds_season.where(ds_season["collection_method"] == 'glider', drop=True)
        elif vers == 'vessel_only':
            ds_season = ds_season.where(ds_season["collection_method"] == 'vessel', drop=True)

        # convert to dataframe and drop nans for the variable you're plotting
        df_season = ds_season.to_dataframe().reset_index()
        df_season = df_season.dropna(subset=[f'{variable}_{stype}'])
        df_season["year"] = df_season["year"].astype(int)
        df_season["month"] = df_season["month"].astype(int)

        years = np.unique(df_season["year"])
        min_year = np.min(years)
        max_year = np.max(years)
        year_list = [y for y in range(min_year, max_year + 1)]

        # plot map of for entire dataset
        fig, ax = cplt.create(extent, **kwargs)
        plt.subplots_adjust(top=.92, bottom=0.08, right=.96, left=0)

        CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                         transform=ccrs.PlateCarree())
        ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

        if clims:
            sct = ax.scatter(df_season['longitude'], df_season['latitude'], c=df_season[f'{variable}_{stype}'], marker='.',
                             vmin=clims[0], vmax=clims[1], s=100, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(),
                             zorder=10)
        else:
            sct = ax.scatter(df_season['longitude'], df_season['latitude'], c=df_season[f'{variable}_{stype}'], marker='.',
                             s=100, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(), zorder=10)

        # Set colorbar height equal to plot height
        divider = make_axes_locatable(ax)
        cax = divider.new_horizontal(size='5%', pad=0.15, axes_class=plt.Axes)
        fig.add_axes(cax)

        if variable == 'omega':
            cbarlab = 'Aragonite Saturation State'
        elif variable == 'pH':
            cbarlab = 'pH'

        # generate colorbar
        cb = plt.colorbar(sct, cax=cax, extend='both')
        cb.set_label(label=cbarlab)

        # add title
        ttl = f'{stype.capitalize()} {variable}: {season} {min_year}-{max_year}'
        if vers == 'glider_only':
            ttl = f'{ttl} (gliders only)'
        elif vers == 'vessel_only':
            ttl = f'{ttl} (vessel only)'
        plt.suptitle(ttl, x=.48, y=.96)

        if addns:  # add NOAA strata to maps
            strata_mapping = cf.return_noaa_polygons()
            for key, values in strata_mapping.items():
                outside_poly = values['poly']
                x, y = outside_poly.exterior.xy
                ax.plot(x, y, color='cyan', lw=2, transform=ccrs.PlateCarree(), zorder=20)
            sfile = os.path.join(savedir, 'with_strata', f'{stype}_{variable}_map_{min_year}-{max_year}-{season.lower()}-{vers}-strata.png')
            os.makedirs(os.path.join(savedir, 'with_strata'), exist_ok=True)
        else:
            sfile = os.path.join(savedir, f'{stype}_{variable}_map_{min_year}-{max_year}-{season.lower()}-{vers}.png')
        plt.savefig(sfile, dpi=200)
        plt.close()

        if vers == 'all':
            # plot maps of summer bottom omega data for each year
            for y in year_list:
                # grab data for the year
                # winter spans two years - grab the December data from the previous year
                if season_code == 'DJF':
                    # Dec from previous year
                    df_season_year1 = df_season[(df_season['year'] == (y - 1)) & (df_season['month'] == 12)]
                    # Jan and Feb for the current year
                    df_season_year2 = df_season[(df_season['year'] == y) & np.logical_or(df_season['month'] == 1, df_season['month'] == 2)]
                    df_season_year = pd.concat([df_season_year1, df_season_year2])
                else:
                    df_season_year = df_season[(df_season['year'] == y) & (df_season['season'] == season_code)]

                fig, ax = cplt.create(extent, **kwargs)
                plt.subplots_adjust(top=.92, bottom=0.08, right=.96, left=0)

                # add bathymetry
                CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                                 transform=ccrs.PlateCarree())
                ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

                if clims:
                    sct = ax.scatter(df_season_year['longitude'], df_season_year['latitude'], c=df_season_year[f'{variable}_{stype}'],
                                     marker='.', vmin=clims[0], vmax=clims[1],
                                     s=100, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(), zorder=10)
                else:
                    sct = ax.scatter(df_season_year['longitude'], df_season_year['latitude'], c=df_season_year[f'{variable}_{stype}'],
                                     marker='.', s=100, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(), zorder=10)

                # Set colorbar height equal to plot height
                divider = make_axes_locatable(ax)
                cax = divider.new_horizontal(size='5%', pad=0.15, axes_class=plt.Axes)
                fig.add_axes(cax)

                # generate colorbar
                cb = plt.colorbar(sct, cax=cax, extend='both')
                cb.set_label(label=cbarlab)

                # add title
                plt.suptitle(f'{stype.capitalize()} {variable}: {season} {y}', x=.48, y=.96)

                savedir_yearly = os.path.join(savedir, f'{stype}_{variable}_years', f'{season.lower()}')
                os.makedirs(savedir_yearly, exist_ok=True)
                sfile = os.path.join(savedir_yearly, f'{stype}_{variable}_map_{y}-{season.lower()}.png')
                
                if addns:  # add NOAA strata to maps
                    strata_mapping = cf.return_noaa_polygons()
                    for key, values in strata_mapping.items():
                        outside_poly = values['poly']
                        x, y = outside_poly.exterior.xy
                        ax.plot(x, y, color='cyan', lw=2, transform=ccrs.PlateCarree(), zorder=20)
                    sfile = os.path.join(savedir_yearly, 'with_strata', f'{stype}_{variable}_map_{min_year}-{max_year}-{season.lower()}-{vers}-strata.png')
                    os.makedirs(os.path.join(savedir_yearly, 'with_strata'), exist_ok=True)

                plt.savefig(sfile, dpi=200)
                plt.close()

                # add to summary
                glider_df = df_season_year[df_season_year['collection_method'] == 'glider']
                glidern = len(glider_df)
                if glidern > 0:
                    glider_date_range = f'{min(glider_df['time']).strftime("%Y%m%d")}-{max(glider_df['time']).strftime("%Y%m%d")}'
                    min_glider_var = np.round(np.nanmin(glider_df[f'{variable}_{stype}']), 2)
                    max_glider_var = np.round(np.nanmax(glider_df[f'{variable}_{stype}']), 2)
                    min_glider_temp = np.round(np.nanmin(glider_df[f'temperature_{stype}']), 2)
                    max_glider_temp = np.round(np.nanmax(glider_df[f'temperature_{stype}']), 2)
                    deployments = np.unique(glider_df['cruise_deployment']).tolist()
                else:
                    glider_date_range = ''
                    min_glider_var = ''
                    max_glider_var = ''
                    min_glider_temp = ''
                    max_glider_temp = ''
                    deployments = ''

                vessel_data = df_season_year[df_season_year['collection_method'] == 'vessel']
                vesseln = len(vessel_data)
                if vesseln > 0:
                    vessel_date_range = f'{min(vessel_data['time']).strftime("%Y%m%d")}-{max(vessel_data['time']).strftime("%Y%m%d")}'
                    min_vessel_var = np.round(np.nanmin(vessel_data[f'{variable}_{stype}']), 2)
                    max_vessel_var = np.round(np.nanmax(vessel_data[f'{variable}_{stype}']), 2)
                    min_vessel_temp = np.round(np.nanmin(vessel_data[f'temperature_{stype}']), 2)
                    max_vessel_temp = np.round(np.nanmax(vessel_data[f'temperature_{stype}']), 2)
                    cruises = np.unique(vessel_data['cruise_deployment']).tolist()
                else:
                    vessel_date_range = ''
                    min_vessel_var = ''
                    max_vessel_var = ''
                    min_vessel_temp = ''
                    max_vessel_temp = ''
                    cruises = ''

                if np.logical_or(glidern > 0, vesseln > 0):
                    rows.append([str(y), season.lower(), glider_date_range, vessel_date_range, str(min_glider_var), str(max_glider_var),
                                 str(min_glider_temp), str(max_glider_temp), glidern, str(min_vessel_var),
                                 str(max_vessel_var), str(min_vessel_temp), str(max_vessel_temp), vesseln, cruises,
                                 deployments])

    summarydf = pd.DataFrame(rows, columns=columns)
    if len(summarydf) > 0:
        summarydf.to_csv(os.path.join(savedir, f'{stype}_{variable}_summary.csv'), index=False)


if __name__ == '__main__':
    stype = 'surface'  # surface bottom
    variable = 'omega'  # pH omega
    version = 'all'  # glider_only vessel_only all
    color_lims = [0.8, 2.2]  # None pH: [7.7, 8.1]  omega: [0.8, 2.2]
    add_noaa_strata = False  # add noaa strata inshore/midshelf/offshore
    save_directory = 'rucool/Saba/NOAA_SOE/data/plots2027'
    main(stype, variable, version, color_lims, add_noaa_strata, save_directory)
