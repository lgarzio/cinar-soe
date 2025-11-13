#!/usr/bin/env python

"""
Author: Lori Garzio on 11/6/2024
Last modified: 11/13/2025
Plot chl/sst from monthly SNPP-VIIRS satellite data downloaded from https://oceancolor.gsfc.nasa.gov/l3/
Calculate seasonal averages from monthly files and plot surface maps.
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


def main(filedir, ab, addns, clims, diff_lims, varname, savedir):
    if ab:
        sdir_avgs = os.path.join(savedir, 'with_boxes', 'seasonal_avgs')
        sdir_anom = os.path.join(savedir, 'with_boxes', 'seasonal_anomalies')
    elif addns:
        sdir_avgs = os.path.join(savedir, 'with_strata', 'seasonal_avgs')
        sdir_anom = os.path.join(savedir, 'with_strata', 'seasonal_anomalies')
    else:
        sdir_avgs = os.path.join(savedir, 'seasonal_avgs')
        sdir_anom = os.path.join(savedir, 'seasonal_anomalies')
    os.makedirs(sdir_avgs, exist_ok=True)
    os.makedirs(sdir_anom, exist_ok=True)

    season_mapping = dict(
        DJF=dict(name='Winter', months=[12, 1, 2]),
        MAM=dict(name='Spring', months=[3, 4, 5]),
        JJA=dict(name='Summer', months=[6, 7, 8]),
        SON=dict(name='Fall', months=[9, 10, 11])
    )

    # regions to subset for average calculations
    nj_box_lons = [-73.9, -72.3, -73.7, -74.8, -74.05, -73.9]
    nj_box_lats = [40.45, 39.9, 38.5, 38.9, 39.65, 40.45]
    gom_box_lons = [-67, -67, -70.3, -70.3, -67]
    gom_box_lats = [44.6, 42.15, 42.15, 43.3, 44.6]

    # # locations of "hot spot" boxes (Gulf of Maine and NJ)
    # gom_box_lons = [-70.4, -69.25, -69.7, -70.85, -70.4]
    # gom_box_lats = [43.36, 43.1, 42.15, 42.4, 43.36]
    # nj_box_lons = [-73.86, -72.75, -73.08, -74.2, -73.86]
    # nj_box_lats = [40.4, 40.14, 39.28, 39.55, 40.4]

    bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    extent = [-78, -65, 35, 45]
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))
    bath_lat = bathy.variables['lat'][:]
    bath_lon = bathy.variables['lon'][:]
    bath_elev = bathy.variables['elevation'][:]

    levels = [-3000, -1000, -100]

    # build a dictionary containing each year-season and which files correspond to that time
    years = np.arange(2012, 2026, 1)
    file_summary_years = dict()

    for y in years:
        file_summary_years[str(y)] = dict(
            DJF=[],
            MAM=[],
            JJA=[],
            SON=[]
        )

    # initialize a dictionary to append all year-season averages
    data_summary_years = dict(
        DJF=dict(),
        MAM=dict(),
        JJA=dict(),
        SON=dict()
    )

    file_list = sorted(glob.glob(os.path.join(filedir, '*.nc')))
    for f in file_list:
        timerange = f.split('.')[1]
        t0 = pd.to_datetime(timerange.split('_')[0])
        t1 = pd.to_datetime(timerange.split('_')[1])

        if t0.month == t1.month:
            for code, sm in season_mapping.items():
                if t0.month in sm['months']:
                    month_code = code
            year = t0.year
            if t0.month == 12:
                year = year + 1
            file_summary_years[str(year)][month_code].append(f)
        else:
            raise(ValueError(f'file contains data from multiple months: {f}'))

    # make a dictionary for overall averages (climatology)
    file_summary = dict(
        DJF=[],
        MAM=[],
        JJA=[],
        SON=[]
    )

    # generate plots for yearly seasonal averages
    for y, items in file_summary_years.items():
        for season_code, files in items.items():
            if len(files) == 3:
                # calculate monthly average
                ds1 = xr.open_dataset(files[0])
                ds2 = xr.open_dataset(files[1])
                ds3 = xr.open_dataset(files[2])

                ds1lat = np.round(ds1.lat.values, 4)
                ds2lat = np.round(ds2.lat.values, 4)
                ds3lat = np.round(ds3.lat.values, 4)
                ds1lon = np.round(ds1.lon.values, 4)
                ds2lon = np.round(ds2.lon.values, 4)
                ds3lon = np.round(ds3.lon.values, 4)

                if np.logical_and(np.sum(ds1lat != ds2lat) == 0, np.sum(ds2lat != ds3lat) == 0):
                    if np.logical_and(np.sum(ds1lon != ds2lon) == 0, np.sum(ds2lon != ds3lon) == 0):
                        # add filenames to the overall summary
                        file_summary[season_code].append(files)

                        data = np.nanmean([ds1[varname].values, ds2[varname].values, ds3[varname].values], axis=0)

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

                        if varname == 'chlor_a':
                            cmap = cmo.cm.algae
                            cmap.set_bad('white')
                            cbar_extend = 'max'
                            cbar_lab = 'Avg Chlorophyll-a (mg m$^{-3}$)'
                            vartitle = 'Average Chlorophyll-a'
                            data = np.ma.masked_less(data, 0.05)
                        elif varname == 'sst_triple':
                            cmap = cmo.cm.thermal
                            cbar_extend = 'both'
                            cbar_lab = 'Average SST (\N{DEGREE SIGN}C)'
                            vartitle = 'Average SST'

                        if clims:
                            bins = clims[1] - clims[0]
                            if varname == 'chlor_a':
                                bins = bins * 2
                            levs = MaxNLocator(nbins=bins).tick_values(clims[0], clims[1])
                            norm = BoundaryNorm(levs, ncolors=cmap.N, clip=True)

                            h = ax.pcolormesh(ds1.lon.values, ds1.lat.values, data, cmap=cmap, norm=norm,
                                              transform=ccrs.PlateCarree(), shading='nearest')
                        else:
                            h = ax.pcolormesh(ds1.lon.values, ds1.lat.values, data, cmap=cmap,
                                              transform=ccrs.PlateCarree(), shading='nearest')

                        # Set colorbar height equal to plot height
                        divider = make_axes_locatable(ax)
                        cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
                        fig.add_axes(cax)

                        # generate colorbar
                        cb = plt.colorbar(h, cax=cax, extend=cbar_extend)
                        cb.set_label(label=cbar_lab)

                        dates = []
                        for f in files:
                            timerange = f.split('.')[1]
                            t0 = pd.to_datetime(timerange.split('_')[0])
                            dates.append(t0)
                        t0str = min(dates).strftime("%b%Y")
                        t1str = max(dates).strftime("%b%Y")

                        # add title
                        ttl = f'{vartitle}: {season} ({t0str}-{t1str})'
                        plt.suptitle(ttl, x=.48, y=.96)

                        if ab:
                            ax.fill(gom_box_lons, gom_box_lats, color='none', edgecolor='magenta', linewidth=2,
                                    transform=ccrs.PlateCarree(), zorder=10)
                            ax.fill(nj_box_lons, nj_box_lats, color='none', edgecolor='magenta', linewidth=2,
                                    transform=ccrs.PlateCarree(), zorder=10)
                        if addns:  # add NOAA strata to maps
                            strata_mapping = cf.return_noaa_polygons()
                            for key, values in strata_mapping.items():
                                outside_poly = values['poly']
                                xx, yy = outside_poly.exterior.xy
                                ax.plot(xx, yy, color='cyan', lw=2, transform=ccrs.PlateCarree(), zorder=20)

                        sfile = os.path.join(sdir_avgs, f'{season.lower()}{str(y)}_{varname}.png')
                        plt.savefig(sfile, dpi=200)
                        plt.close()

    # calculate seasonal climatologies (averages) for all years
    for season_code, nested_files in file_summary.items():
        files = list(itertools.chain(*nested_files))
        dates = []
        for i, f in enumerate(files):
            timerange = f.split('.')[1]
            t0 = pd.to_datetime(timerange.split('_')[0])
            dates.append(t0)

            ds = xr.open_dataset(f)
            if i == 0:
                lat = ds.lat.values
                lon = ds.lon.values
                data = np.empty((len(files), len(lat), len(lon)))  # set up the empty 3D array to append data
                data[:] = np.nan
            data[i] = ds[varname].values  # add data to the correct spot in the matrix

        # calculate average for all years (seasonal climatology)
        climatology = np.nanmean(data, axis=0)
        t0str = min(dates).strftime("%b%Y")
        t1str = max(dates).strftime("%b%Y")

        season = season_mapping[season_code]['name']

        # plot maps of overall seasonal averages
        fig, ax = plt.subplots(figsize=(9, 8), subplot_kw=dict(projection=ccrs.Mercator()))
        plt.subplots_adjust(top=.92, bottom=0.08, right=.96, left=0)

        # define bathymetry levels and data
        CS = plt.contour(bath_lon, bath_lat, bath_elev, levels, linewidths=.75, alpha=.5, colors='k',
                         transform=ccrs.PlateCarree())
        ax.clabel(CS, levels, inline=True, fontsize=7, fmt='%d')

        cf.add_map_features(ax, extent)

        if varname == 'chlor_a':
            cmap = cmo.cm.algae
            cmap.set_bad('white')
            cbar_extend = 'max'
            cbar_lab = 'Avg Chlorophyll-a (mg m$^{-3}$)'
            vartitle = 'Average Chlorophyll-a'
            data = np.ma.masked_less(data, 0.05)
        elif varname == 'sst_triple':
            cmap = cmo.cm.thermal
            cbar_extend = 'both'
            cbar_lab = 'Average SST (\N{DEGREE SIGN}C)'
            vartitle = 'Average SST'

        if clims:
            bins = clims[1] - clims[0]
            if varname == 'chlor_a':
                bins = bins * 2
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
            # ax.fill(gom_box_lons, gom_box_lats, color='none', edgecolor='cyan', linewidth=2.5,
            #         transform=ccrs.PlateCarree(), zorder=10)
            ax.fill(nj_box_lons, nj_box_lats, color='none', edgecolor='magenta', linewidth=2,
                    transform=ccrs.PlateCarree(), zorder=10)
        if addns:  # add NOAA strata to maps
            strata_mapping = cf.return_noaa_polygons()
            for key, values in strata_mapping.items():
                outside_poly = values['poly']
                xx, yy = outside_poly.exterior.xy
                ax.plot(xx, yy, color='cyan', lw=2, transform=ccrs.PlateCarree(), zorder=20)

        sfile = os.path.join(sdir_avgs, f'{season.lower()}-{varname}-climatology.png')
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

            if varname == 'chlor_a':
                cbar_lab = 'Chlorophyll-a Anomaly (mg m$^{-3}$)'
                vartitle = 'Chl-a Anomaly'
                masked_diff = np.ma.masked_inside(diff, -.1, .1)
            elif varname == 'sst_triple':
                cbar_lab = 'SST Anomaly (\N{DEGREE SIGN}C)'
                vartitle = 'SST Anomaly'
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
                # ax.fill(gom_box_lons, gom_box_lats, color='none', edgecolor='magenta', linewidth=2,
                #         transform=ccrs.PlateCarree(), zorder=10)
                ax.fill(nj_box_lons, nj_box_lats, color='none', edgecolor='magenta', linewidth=2,
                        transform=ccrs.PlateCarree(), zorder=10)
            if addns:  # add NOAA strata to maps
                strata_mapping = cf.return_noaa_polygons()
                for key, values in strata_mapping.items():
                    outside_poly = values['poly']
                    xx, yy = outside_poly.exterior.xy
                    ax.plot(xx, yy, color='cyan', lw=2, transform=ccrs.PlateCarree(), zorder=20)

            sfile = os.path.join(sdir_anom, f'{season.lower()}{year}-{varname}-anomaly.png')
            plt.savefig(sfile, dpi=200)
            plt.close()


if __name__ == '__main__':
    filedir = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/nasa_data/chlor_a/opendap_monthly'
    #filedir = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/nasa_data/sst_triple/opendap'
    add_boxes = False
    add_noaa_strata = False  # add noaa strata inshore/midshelf/offshore
    color_lims = [0, 10]  # None sst: [0, 28] chla: [0, 10]
    diff_lims = [-4, 4, 20]  # [min, max, bins] for difference plots  sst: [-5, 5, 40] chla: [-4, 4, 20]
    varname = 'chlor_a'  # sst_triple chlor_a
    save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/plots2026/chlor_a'
    # save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/plots2025/sst'
    main(filedir, add_boxes, add_noaa_strata, color_lims, diff_lims, varname, save_directory)
