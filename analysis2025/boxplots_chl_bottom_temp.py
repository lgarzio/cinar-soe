#!/usr/bin/env python

"""
Author: Lori Garzio on 12/13/2024
Last modified: 12/13/2024
Generate boxplots for satellite-derived chla and GLORYS bottom temperature for summers of 2022-2024 for NYB inshore,
midshelf, and offshore regions defined by modified NOAA bottom trawl survey strata downloaded from
https://github.com/NOAA-EDAB/FisheryConditionLinks/tree/master/NES_BOTTOM_TRAWL_STRATA.
The box limits extend from the lower to upper quartiles (25%, 75%), with a line at the median and a diamond symbol at
the mean. The whiskers extend from the box by 1.5x the inter-quartile range (IQR). Circles indicate outliers.
Notch indicates 95% CI around the median.
"""

import os
import glob
import pandas as pd
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console
plt.rcParams.update({'font.size': 15})


def set_box_colors(bp, colors):
    for key in ['boxes', 'medians', 'fliers', 'means']:
        for patch, color in zip(bp[key], colors):
            patch.set_color(color)
            if key == 'boxes':
                patch.set_facecolor('none')
            elif key == 'means':
                patch.set_markerfacecolor(color)
                patch.set_markeredgecolor(color)
            elif key == 'fliers':
                patch.set_markeredgecolor(color)

    wc_colors = [x for pair in zip(colors, colors) for x in pair]
    for key in ['whiskers', 'caps']:
        for patch, color in zip(bp[key], wc_colors):
            patch.set_color(color)


def main(varname, savedir):
    strata_mapping = cf.return_noaa_polygons()

    datadict = dict(
        inshore=dict(),
        midshelf=dict(),
        offshore=dict()
    )

    colors = dict(
        y2022='tab:blue',
        y2023='#d7369eff',
        y2024='k'
    )

    if varname == 'chlor_a':
        ylab = 'Monthly Chlorophyll-a (mg m$^{-3}$)'
        ymult = 0

        filedir = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/nasa_data/surface_chl/requested_files'
        file_list = sorted(glob.glob(os.path.join(filedir, '*.nc')))
        for f in file_list:
            timerange = f.split('.')[1]
            t0 = pd.to_datetime(timerange.split('_')[0])
            fyear = t0.year

            if fyear in [2022, 2023, 2024]:
                if t0.month in [6, 7, 8]:
                    ds = xr.open_dataset(f)

                    # select data inside NOAA bottom trawl survey polygons for inshore/midshelf/offshore
                    for code, sm in strata_mapping.items():
                        outside_poly = sm['poly']
                        x, y = outside_poly.exterior.xy

                        # select data inside the polygon
                        long, latg = np.meshgrid(ds.lon, ds.lat)
                        long, latg = long.flatten(), latg.flatten()
                        points = np.vstack((long, latg)).T

                        polygon_bounds = list(zip(x, y))
                        bounds = Path(polygon_bounds)
                        grid = ~bounds.contains_points(points)
                        grid_mask = grid.reshape(len(ds.lat), len(ds.lon))
                        data = ds.chlor_a.values.copy()
                        data[grid_mask] = np.nan

                        data = data.flatten()
                        data = data[~np.isnan(data)]

                        # add to dictionary
                        try:
                            datadict[code][str(fyear)] = np.append(datadict[code][str(fyear)], data)
                        except KeyError:
                            datadict[code][str(fyear)] = data

    elif varname == 'bottom_temp':
        ylab = 'Monthly Bottom Temperature (\N{DEGREE SIGN}C)'
        ymult = 1.01

        # files are grouped by year-season so one file will have 3 months of data
        filedir = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/cmems_data/glorys'
        file_list = sorted(glob.glob(os.path.join(filedir, '*.nc')))

        for f in file_list:
            start_str = f.split('-')[-1].split('_')[0]
            t0 = pd.to_datetime(start_str)
            fyear = t0.year

            if fyear in [2022, 2023, 2024]:
                if t0.month == 6:
                    ds = xr.open_dataset(f)

                    for code, sm in strata_mapping.items():
                        outside_poly = sm['poly']
                        x, y = outside_poly.exterior.xy

                        # select data inside the polygon
                        long, latg = np.meshgrid(ds.longitude, ds.latitude)
                        long, latg = long.flatten(), latg.flatten()
                        points = np.vstack((long, latg)).T

                        polygon_bounds = list(zip(x, y))
                        bounds = Path(polygon_bounds)
                        grid = ~bounds.contains_points(points)
                        grid_mask = grid.reshape(len(ds.latitude), len(ds.longitude))
                        for t in ds.time:
                            data = ds.sel(time=t)['bottomT'].values.copy()
                            data[grid_mask] = np.nan

                            data = data.flatten()
                            data = data[~np.isnan(data)]

                            # add to dictionary
                            try:
                                datadict[code][str(fyear)] = np.append(datadict[code][str(fyear)], data)
                            except KeyError:
                                datadict[code][str(fyear)] = data

    # for boxplots
    inshore = [list(datadict['inshore']['2022']),
               list(datadict['inshore']['2023']),
               list(datadict['inshore']['2024'])]
    midshelf = [list(datadict['midshelf']['2022']),
                list(datadict['midshelf']['2023']),
                list(datadict['midshelf']['2024'])]
    offshore = [list(datadict['offshore']['2022']),
                list(datadict['offshore']['2023']),
                list(datadict['offshore']['2024'])]

    fig, ax = plt.subplots(figsize=(9, 8))
    plt.subplots_adjust(top=.94, bottom=0.08, right=.96, left=0.1)

    # customize the boxplot elements
    meanpointprops = dict(marker='D')
    box_colors = [colors['y2022'], colors['y2023'], colors['y2024']]

    bp_in = ax.boxplot(inshore, positions=[1, 2, 3], widths=0.6, patch_artist=True, showmeans=True,
                       notch=True, meanprops=meanpointprops, sym='.')
    bp_mid = ax.boxplot(midshelf, positions=[5, 6, 7], widths=0.6, patch_artist=True, showmeans=True,
                        notch=True, meanprops=meanpointprops, sym='.')
    bp_off = ax.boxplot(offshore, positions=[9, 10, 11], widths=0.6, patch_artist=True, showmeans=True,
                        notch=True, meanprops=meanpointprops, sym='.')

    # set box colors
    set_box_colors(bp_in, box_colors)
    set_box_colors(bp_mid, box_colors)
    set_box_colors(bp_off, box_colors)

    # draw temporary lines and use them to create a legend
    plt.plot([], c=colors['y2022'], label='2022')
    plt.plot([], c=colors['y2023'], label='2023')
    plt.plot([], c=colors['y2024'], label='2024')
    plt.legend(fontsize=12)

    # set axes labels
    ax.set_xticks([2, 6, 10])
    ax.set_xticklabels(['Inshore', 'Midshelf', 'Offshore'])
    ax.set_ylabel(ylab)

    ylims = ax.get_ylim()
    ax.set_ylim(ylims)
    if varname == 'chlor_a':
        ax.set_ylim(-.5, 14.5)  # chla
    ax.set_xlim(0, 12)
    ax.vlines(4, ylims[0], ylims[1], colors='k', lw=.8)
    ax.vlines(8, ylims[0], ylims[1], colors='k', lw=.8)

    # add sample sizes for each box
    pos = [1, 2, 3]
    for ii, jj in enumerate(inshore):
        ss = f'{len(jj)}'
        ax.text(pos[ii], ylims[0] * ymult, ss, size=9, horizontalalignment='center')

    pos = [5, 6, 7]
    for ii, jj in enumerate(midshelf):
        ss = f'{len(jj)}'
        ax.text(pos[ii], ylims[0] * ymult, ss, size=9, horizontalalignment='center')

    pos = [9, 10, 11]
    for ii, jj in enumerate(offshore):
        ss = f'{len(jj)}'
        ax.text(pos[ii], ylims[0] * ymult, ss, size=9, horizontalalignment='center')

    sfile = os.path.join(savedir, f'{varname}_boxplot_summer2022-2024.png')
    plt.savefig(sfile, dpi=300)
    plt.close()


if __name__ == '__main__':
    variable = 'bottom_temp'  # chlor_a bottom_temp
    save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/plots2025/boxplots'
    main(variable, save_directory)
