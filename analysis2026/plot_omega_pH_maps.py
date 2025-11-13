#!/usr/bin/env python

"""
Author: Lori Garzio on 10/18/2021
Last modified: 11/13/2025
Plot bottom/surface omega/pH maps using CODAP-NA, EcoMon, and glider datasets.
CODAP-NA dataset documented here: https://essd.copernicus.org/articles/13/2777/2021/
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


def main(cruise_file, glider_file, ab, addns, stype, variable, vers, clims, savedir):
    savedir = os.path.join(savedir, f'{stype}_{variable}')
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

    bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    extent = [-78, -65, 35, 45]
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))

    cdata = xr.open_dataset(cruise_file)
    gdata = xr.open_dataset(glider_file)

    # add month and season
    cdata['month'] = cdata['time.month']
    cdata['season'] = cdata['time.season']
    gdata['month'] = gdata['time.month']
    gdata['season'] = gdata['time.season']

    # convert to dataframe and merge
    cdf = cdata.to_dataframe()
    gdf = gdata.to_dataframe()
    cdf['deployment'] = ''
    gdf['cruise'] = ''
    cdf['sample_type'] = 'vessel'
    gdf['sample_type'] = 'glider'
    cdf = cdf[['cruise', 'deployment', 'lat', 'lon', 'depth', f'pH_{stype}', f'omega_{stype}', f'temperature_{stype}', 'month', 'season', 'sample_type']]
    gdf = gdf[['cruise', 'deployment', 'lat', 'lon', 'depth', f'pH_{stype}', f'omega_{stype}', f'temperature_{stype}', 'month', 'season', 'sample_type']]
    df = pd.concat([cdf, gdf])

    # regions to subset for average calculations
    nj_box_lons = [-73.9, -71.7, -73.4, -74.8, -74.05, -73.9]
    nj_box_lats = [40.45, 39.7, 38.3, 38.9, 39.65, 40.45]
    # nj_box_lons2 = [-72, -71.6, -73.7, -74.95, -74.05, -73.9, -72]  # larger NYB region
    # nj_box_lats2 = [40.9, 39.7, 38, 38.7, 39.65, 40.45, 40.9]  # larger NYB region
    # gom_box_lons = [-67, -67, -70.5, -70.5, -70, -69, -67]
    # gom_box_lats = [44.6, 42.2, 42.2, 43.05, 43.55, 44, 44.6]

    for season_code in np.unique(gdata.season):
        print(season_code)
        season = season_mapping[season_code]
        df_season = df[df.season == season_code]

        if vers == 'glider_only':
            df_season = df_season[df.sample_type == 'glider']
        elif vers == 'vessel_only':
            df_season = df_season[df.sample_type == 'vessel']

        years = np.unique(pd.to_datetime(df_season.index).year)
        min_year = np.min(years)
        max_year = np.max(years)
        year_list = [y for y in range(min_year, max_year + 1)]

        # plot map of for entire dataset
        fig, ax = plt.subplots(figsize=(9, 8), subplot_kw=dict(projection=ccrs.Mercator()))
        plt.subplots_adjust(top=.92, bottom=0.08, right=.96, left=0)

        # define bathymetry levels and data
        bath_lat = bathy['lat'].values
        bath_lon = bathy['lon'].values
        bath_elev = bathy['elevation'].values

        levels = [-3000, -1000, -100]
        CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                         transform=ccrs.PlateCarree())
        ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

        #cplt.add_bathymetry(ax, bath_lon, bath_lat, bath_elev, method='topo_log', zorder=1)

        cf.add_map_features(ax, extent, oceancolor=cfeature.COLORS['water'])

        if clims:
            sct = ax.scatter(df_season['lon'], df_season['lat'], c=df_season[f'{variable}_{stype}'], marker='.',
                             vmin=clims[0], vmax=clims[1], s=100, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(),
                             zorder=10)
        else:
            sct = ax.scatter(df_season['lon'], df_season['lat'], c=df_season[f'{variable}_{stype}'], marker='.',
                             s=100, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(), zorder=10)

        # Set colorbar height equal to plot height
        divider = make_axes_locatable(ax)
        cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
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

        if ab:
            # ax.fill(gom_box_lons, gom_box_lats, color='none', edgecolor='cyan', linewidth=3,
            #         transform=ccrs.PlateCarree(), zorder=10)
            # ax.fill(nj_box_lons2, nj_box_lats2, color='none', edgecolor='magenta', linewidth=3,
            #         transform=ccrs.PlateCarree(), zorder=10)
            ax.fill(nj_box_lons, nj_box_lats, color='none', edgecolor='cyan', linewidth=3,
                    transform=ccrs.PlateCarree(), zorder=10)

            sfile = os.path.join(savedir, 'with_boxes', f'{stype}_{variable}_map_{min_year}-{max_year}-{season.lower()}-{vers}-boxes.png')
            os.makedirs(os.path.join(savedir, 'with_boxes'), exist_ok=True)
        elif addns:  # add NOAA strata to maps
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
                    df_season_year1 = df_season[(df_season.index.year == (y - 1)) & (df_season['month'] == 12)]
                    # Jan and Feb for the current year
                    df_season_year2 = df_season[(df_season.index.year == y) & np.logical_or(df_season['month'] == 1, df_season['month'] == 2)]
                    df_season_year = pd.concat([df_season_year1, df_season_year2])
                else:
                    df_season_year = df_season[(df_season.index.year == y) & (df_season['season'] == season_code)]

                fig, ax = plt.subplots(figsize=(9, 8), subplot_kw=dict(projection=ccrs.Mercator()))
                plt.subplots_adjust(top=.92, bottom=0.08, right=.96, left=0)

                # add bathymetry
                CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                                 transform=ccrs.PlateCarree())
                ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

                #cplt.add_bathymetry(ax, bath_lon, bath_lat, bath_elev, method='topo_log', zorder=1)

                cf.add_map_features(ax, extent, oceancolor=cfeature.COLORS['water'])
                # cf.add_map_features(ax, extent)

                if clims:
                    sct = ax.scatter(df_season_year['lon'], df_season_year['lat'], c=df_season_year[f'{variable}_{stype}'],
                                     marker='.', vmin=clims[0], vmax=clims[1],
                                     s=100, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(), zorder=10)
                else:
                    sct = ax.scatter(df_season_year['lon'], df_season_year['lat'], c=df_season_year[f'{variable}_{stype}'],
                                     marker='.', s=100, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(), zorder=10)

                # Set colorbar height equal to plot height
                divider = make_axes_locatable(ax)
                cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
                fig.add_axes(cax)

                # generate colorbar
                cb = plt.colorbar(sct, cax=cax, extend='both')
                cb.set_label(label=cbarlab)

                # add title
                plt.suptitle(f'{stype.capitalize()} {variable}: {season} {y}', x=.48, y=.96)

                if ab:
                    # ax.fill(gom_box_lons, gom_box_lats, color='none', edgecolor='cyan', linewidth=3,
                    #         transform=ccrs.PlateCarree(), zorder=10)
                    # ax.fill(nj_box_lons, nj_box_lats, color='none', edgecolor='cyan', linewidth=3,
                    #         transform=ccrs.PlateCarree(), zorder=10)
                    ax.fill(nj_box_lons, nj_box_lats, color='none', edgecolor='cyan', linewidth=3,
                            transform=ccrs.PlateCarree(), zorder=10)
                elif addns:  # add NOAA strata to maps
                    strata_mapping = cf.return_noaa_polygons()
                    for key, values in strata_mapping.items():
                        outside_poly = values['poly']
                        x, y = outside_poly.exterior.xy
                        ax.plot(x, y, color='cyan', lw=2, transform=ccrs.PlateCarree(), zorder=20)

                if len(df_season_year) > 0:
                    versions = ['b', 'c', 'd', 'e', 'f', 'g']
                else:
                    versions = ['b', 'c', 'd']

                savedir_yearly = os.path.join(savedir, f'{stype}_{variable}_years', f'{season.lower()}')
                os.makedirs(savedir_yearly, exist_ok=True)
                for v in versions:
                    if ab:
                        sfile = os.path.join(savedir_yearly, 'with_boxes', f'{stype}_{variable}_map_{y}_{season.lower()}-{v}-boxes.png')
                        os.makedirs(os.path.join(savedir_yearly, 'with_boxes'), exist_ok=True)
                    elif addns:
                        sfile = os.path.join(savedir_yearly, 'with_strata',
                                             f'{stype}_{variable}_map_{min_year}-{max_year}-{season.lower()}-{vers}-strata.png')
                        os.makedirs(os.path.join(savedir_yearly, 'with_strata'), exist_ok=True)
                    else:
                        sfile = os.path.join(savedir_yearly, f'{stype}_{variable}_map_{y}-{season.lower()}-{v}.png')
                    plt.savefig(sfile, dpi=200)
                plt.close()

                # add to summary
                glider_df = df_season_year[df_season_year['sample_type'] == 'glider']
                glidern = len(glider_df)
                if glidern > 0:
                    glider_date_range = f'{min(glider_df.index).strftime("%Y%m%d")}-{max(glider_df.index).strftime("%Y%m%d")}'
                    min_glider_var = np.round(np.nanmin(glider_df[f'{variable}_{stype}']), 2)
                    max_glider_var = np.round(np.nanmax(glider_df[f'{variable}_{stype}']), 2)
                    min_glider_temp = np.round(np.nanmin(glider_df[f'temperature_{stype}']), 2)
                    max_glider_temp = np.round(np.nanmax(glider_df[f'temperature_{stype}']), 2)
                    deployments = np.unique(glider_df['deployment']).tolist()
                else:
                    glider_date_range = ''
                    min_glider_var = ''
                    max_glider_var = ''
                    min_glider_temp = ''
                    max_glider_temp = ''
                    deployments = ''

                vessel_data = df_season_year[df_season_year['sample_type'] == 'vessel']
                vesseln = len(vessel_data)
                if vesseln > 0:
                    vessel_date_range = f'{min(vessel_data.index).strftime("%Y%m%d")}-{max(vessel_data.index).strftime("%Y%m%d")}'
                    min_vessel_var = np.round(np.nanmin(vessel_data[f'{variable}_{stype}']), 2)
                    max_vessel_var = np.round(np.nanmax(vessel_data[f'{variable}_{stype}']), 2)
                    min_vessel_temp = np.round(np.nanmin(vessel_data[f'temperature_{stype}']), 2)
                    max_vessel_temp = np.round(np.nanmax(vessel_data[f'temperature_{stype}']), 2)
                    cruises = np.unique(vessel_data['cruise']).tolist()
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
    cruise_data = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/output_nc/vessel_based_bottom_OA_data_2007_2023.nc'  # bottom
    glider_data = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/output_nc/glider_based_bottom_OA_data_2019_2025.nc'  # bottom
    #cruise_data = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/output_nc/vessel_based_surface_OA_data_2004_2023.nc'  # surface
    #glider_data = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/output_nc/glider_based_surface_OA_data_2019_2025.nc'  # surface
    add_boxes = False
    add_noaa_strata = False  # add noaa strata inshore/midshelf/offshore
    stype = 'bottom'  # surface bottom
    variable = 'omega'  # pH omega
    version = 'all'  # glider_only vessel_only all
    color_lims = [0.8, 2.2]  # None pH: [7.7, 8.1]  omega: [0.8, 2.2]
    save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/plots2026'
    main(cruise_data, glider_data, add_boxes, add_noaa_strata, stype, variable, version, color_lims, save_directory)
