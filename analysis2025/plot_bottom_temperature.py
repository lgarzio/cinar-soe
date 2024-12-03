#!/usr/bin/env python

"""
Author: Lori Garzio on 11/14/2024
Last modified: 11/14/2024
Plot bottom temperature data from GLORYS12V1 (Global Ocean Physics Reanalysis)
https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description
Monthly averaged files were downloaded (3 months per file = 1 season), and seasonal averages were calculated and plotted.
Also calculate seasonal climatologies and plot seasonal-yearly anomalies
"""

import os
import glob
import itertools
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import BoundaryNorm
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import cmocean as cmo
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console
plt.rcParams.update({'font.size': 15})


def main(filedir, ab, stype, clims, diff_lims, savedir):
    sdir_avgs = os.path.join(savedir, 'seasonal_avgs')
    sdir_anom = os.path.join(savedir, 'seasonal_anomalies')

    season_mapping = dict(
        DJF=dict(name='Winter', months=[12, 1, 2]),
        MAM=dict(name='Spring', months=[3, 4, 5]),
        JJA=dict(name='Summer', months=[6, 7, 8]),
        SON=dict(name='Fall', months=[9, 10, 11])
    )

    # locations of "hot spot" boxes (Gulf of Maine and NJ)
    gom_box_lons = [-70.4, -69.25, -69.7, -70.85, -70.4]
    gom_box_lats = [43.36, 43.1, 42.15, 42.4, 43.36]
    nj_box_lons = [-73.86, -72.75, -73.08, -74.2, -73.86]
    nj_box_lats = [40.4, 40.14, 39.28, 39.55, 40.4]

    bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    extent = [-78, -65, 35, 45]
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))
    bath_lat = bathy.variables['lat'][:]
    bath_lon = bathy.variables['lon'][:]
    bath_elev = bathy.variables['elevation'][:]

    levels = [-3000, -1000, -100]

    # initialize a dictionary to append all year-season averages
    data_summary_years = dict(
        DJF=dict(),
        MAM=dict(),
        JJA=dict(),
        SON=dict()
    )

    # make a dictionary for overall averages (climatology)
    file_summary = dict(
        DJF=[],
        MAM=[],
        JJA=[],
        SON=[]
    )

    # generate plots for yearly seasonal averages
    # each file contains a season of data that needs to be averaged
    file_list = sorted(glob.glob(os.path.join(filedir, '*.nc')))
    for f in file_list:
        # calculate monthly average
        ds = xr.open_dataset(f)

        y = int(np.max(np.unique(ds['time.year'].values)))
        season_code = np.unique(ds['time.season'].values)
        if len(season_code) > 1:
            raise ValueError(f'file {f} contains >1 season')
        season_code = season_code[0]

        # add filenames to the overall summary
        file_summary[season_code].append(f)

        data = np.nanmean(ds['bottomT'].values, axis=0)

        # add data to the summary dictionary for later comparison to overall climatology
        data_summary_years[season_code][y] = data

        season = season_mapping[season_code]['name']

        # plot maps of avg seasonal chl for each year
        fig, ax = plt.subplots(figsize=(9, 8), subplot_kw=dict(projection=ccrs.Mercator()))
        plt.subplots_adjust(top=.92, bottom=0.08, right=.96, left=0)

        # define bathymetry levels and data
        CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                         transform=ccrs.PlateCarree())
        ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

        cf.add_map_features(ax, extent)

        cmap = cmo.cm.thermal
        cbar_extend = 'both'
        cbar_lab = 'Average Bottom Temperature (\N{DEGREE SIGN}C)'
        vartitle = 'Average Bottom Temperature'

        if clims:
            bins = clims[1] - clims[0]
            levs = MaxNLocator(nbins=bins).tick_values(clims[0], clims[1])
            norm = BoundaryNorm(levs, ncolors=cmap.N, clip=True)

            h = ax.pcolormesh(ds.longitude.values, ds.latitude.values, data, cmap=cmap, norm=norm,
                              transform=ccrs.PlateCarree(), shading='nearest')
        else:
            h = ax.pcolormesh(ds.longitude.values, ds.latitude.values, data, cmap=cmap,
                              transform=ccrs.PlateCarree(), shading='nearest')

        # Set colorbar height equal to plot height
        divider = make_axes_locatable(ax)
        cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
        fig.add_axes(cax)

        # generate colorbar
        cb = plt.colorbar(h, cax=cax, extend=cbar_extend)
        cb.set_label(label=cbar_lab)

        t0str = min(pd.to_datetime(ds.time.values)).strftime("%b%Y")
        t1str = max(pd.to_datetime(ds.time.values)).strftime("%b%Y")

        # add title
        ttl = f'{vartitle}: {season} ({t0str}-{t1str})'
        plt.suptitle(ttl, x=.48, y=.96)

        if ab:
            ax.fill(gom_box_lons, gom_box_lats, color='none', edgecolor='cyan', linewidth=2.5,
                    transform=ccrs.PlateCarree(), zorder=10)
            ax.fill(nj_box_lons, nj_box_lats, color='none', edgecolor='cyan', linewidth=2.5,
                    transform=ccrs.PlateCarree(), zorder=10)

            sfile = os.path.join(sdir_avgs, f'{season.lower()}{str(y)}_bottom_temp-boxes.png')
        else:
            sfile = os.path.join(sdir_avgs, f'{season.lower()}{str(y)}_bottom_temp.png')
        plt.savefig(sfile, dpi=200)
        plt.close()

    # calculate seasonal climatologies (averages) for all years
    for season_code, files in file_summary.items():
        dates = np.array([], dtype='datetime64[ns]')
        for i, f in enumerate(files):
            ds = xr.open_dataset(f)
            dates = np.append(dates, ds.time.values)
            if i == 0:
                lat = ds.latitude.values
                lon = ds.longitude.values
                data = ds['bottomT'].values
            else:
                data = np.concatenate([data, ds['bottomT'].values])

        # calculate average for all years (seasonal climatology)
        climatology = np.nanmean(data, axis=0)
        t0str = min(pd.to_datetime(dates)).strftime("%b%Y")
        t1str = max(pd.to_datetime(dates)).strftime("%b%Y")

        season = season_mapping[season_code]['name']

        # plot maps of overall seasonal averages
        fig, ax = plt.subplots(figsize=(9, 8), subplot_kw=dict(projection=ccrs.Mercator()))
        plt.subplots_adjust(top=.92, bottom=0.08, right=.96, left=0)

        # define bathymetry levels and data
        CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                         transform=ccrs.PlateCarree())
        ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

        cf.add_map_features(ax, extent)

        cmap = cmo.cm.thermal
        cbar_extend = 'both'
        cbar_lab = 'Average Bottom Temperature (\N{DEGREE SIGN}C)'
        vartitle = 'Average Bottom Temperature'

        if clims:
            bins = clims[1] - clims[0]
            levs = MaxNLocator(nbins=bins).tick_values(clims[0], clims[1])
            norm = BoundaryNorm(levs, ncolors=cmap.N, clip=True)

            h = ax.pcolormesh(lon, lat, climatology, cmap=cmap, norm=norm,
                              transform=ccrs.PlateCarree(), shading='nearest')
        else:
            h = ax.pcolormesh(lon, lat, climatology, cmap=cmap,
                              transform=ccrs.PlateCarree(), shading='nearest')

        # Set colorbar height equal to plot height
        divider = make_axes_locatable(ax)
        cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
        fig.add_axes(cax)

        # generate colorbar
        cb = plt.colorbar(h, cax=cax, extend=cbar_extend)
        cb.set_label(label=cbar_lab)

        # add title
        ttl = f'{vartitle}: {season} ({t0str}-{t1str})'
        plt.suptitle(ttl, x=.48, y=.96)

        if ab:
            ax.fill(gom_box_lons, gom_box_lats, color='none', edgecolor='cyan', linewidth=2.5,
                    transform=ccrs.PlateCarree(), zorder=10)
            ax.fill(nj_box_lons, nj_box_lats, color='none', edgecolor='cyan', linewidth=2.5,
                    transform=ccrs.PlateCarree(), zorder=10)

            sfile = os.path.join(sdir_avgs, f'{season.lower()}-bottom_temp-climatology-boxes.png')
        else:
            sfile = os.path.join(sdir_avgs, f'{season.lower()}-bottom_temp-climatology.png')
        plt.savefig(sfile, dpi=200)
        plt.close()

        # go through each season-year again and compare to the climatology (anomalies)
        for year, yearly_data in data_summary_years[season_code].items():
            diff = yearly_data - climatology

            cmap = plt.get_cmap('RdBu_r')
            cmap.set_bad('white')

            # plot maps of seasonal-yearly anomalies
            fig, ax = plt.subplots(figsize=(9, 8), subplot_kw=dict(projection=ccrs.Mercator()))
            plt.subplots_adjust(top=.92, bottom=0.08, right=.96, left=0)

            # define bathymetry levels and data
            CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                             transform=ccrs.PlateCarree())
            ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

            cf.add_map_features(ax, extent)

            cbar_lab = 'Bottom Temperature Anomaly (\N{DEGREE SIGN}C)'
            vartitle = 'Bottom Temperature Anomaly'
            masked_diff = np.ma.masked_inside(diff, -.1, .1)

            levs = MaxNLocator(nbins=diff_lims[2]).tick_values(diff_lims[0], diff_lims[1])
            norm = BoundaryNorm(levs, ncolors=cmap.N, clip=True)

            h = ax.pcolormesh(lon, lat, masked_diff, cmap=cmap, norm=norm,
                              transform=ccrs.PlateCarree(), shading='nearest')

            # Set colorbar height equal to plot height
            divider = make_axes_locatable(ax)
            cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
            fig.add_axes(cax)

            # generate colorbar
            cb = plt.colorbar(h, cax=cax, extend='both')
            cb.set_label(label=cbar_lab)

            # add title
            ttl = f'{vartitle}: {season} {year} (vs {t0str}-{t1str})'
            plt.suptitle(ttl, x=.48, y=.96)

            if ab:
                ax.fill(gom_box_lons, gom_box_lats, color='none', edgecolor='cyan', linewidth=2.5,
                        transform=ccrs.PlateCarree(), zorder=10)
                ax.fill(nj_box_lons, nj_box_lats, color='none', edgecolor='cyan', linewidth=2.5,
                        transform=ccrs.PlateCarree(), zorder=10)

                sfile = os.path.join(sdir_anom, f'{season.lower()}{year}-bottom_temp-anomaly-boxes.png')
            else:
                sfile = os.path.join(sdir_anom, f'{season.lower()}{year}-bottom_temp-anomaly.png')
            plt.savefig(sfile, dpi=200)
            plt.close()


if __name__ == '__main__':
    filedir = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/cmems_data/glorys'
    add_boxes = False
    stype = 'bottom'
    color_lims = [2, 22]  # None bottom temp: [2, 22]
    diff_lims = [-5, 5, 40]  # [min, max, bins] for difference plots
    save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/plots2025/bottom_temp'
    main(filedir, add_boxes, stype, color_lims, diff_lims, save_directory)
