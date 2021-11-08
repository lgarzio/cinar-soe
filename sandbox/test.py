#!/usr/bin/env python

"""
Author: Lori Garzio on 10/18/2021
Last modified: 11/8/2021
CODAP-NA, glider and EcoMon datasets for OA analysis.
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

    # get CODAP data
    ds = xr.open_dataset(f)
    years = ds.Year_UTC.values.astype(int)
    months = ds.Month_UTC.values.astype(int)
    days = ds.Day_UTC.values.astype(int)
    times = np.array([])
    for i, year in enumerate(years):
        times = np.append(times, dt.datetime(year, months[i], days[i]).strftime('%Y-%m-%dT00:00:00'))

    times = times.astype('datetime64[s]')

    idx = []
    data = dict(Month_UTC=np.array([]),
                Year_UTC=np.array([]),
                Cruise_ID=np.array([]),
                Profile_number=np.array([]),
                Latitude=np.array([]),
                Longitude=np.array([]),
                Depth=np.array([]),
                Depth_bottom=np.array([]),
                CTDTEMP_ITS90=np.array([]),
                CTDTEMP_flag=np.array([]),
                recommended_Salinity_PSS78=np.array([]),
                recommended_Salinity_flag=np.array([]),
                pH_TS_insitu_calculated=np.array([]),
                pH_TS_insitu_measured=np.array([]),
                Aragonite=np.array([]))
    for i, lon in enumerate(ds.Longitude.values):
        if Polygon(list(zip(lon_bounds, lat_bounds))).contains(Point(lon, ds.Latitude.values[i])):
            idx.append(i)
            for key in data.keys():
                if key == 'Cruise_ID':
                    cid = ds.Cruise_ID.values[:, i]
                    cid = [x.decode('UTF-8') for x in cid]
                    data[key] = np.append(data[key], ''.join(cid).strip())
                else:
                    data[key] = np.append(data[key], ds[key].values[i])

    # select where bottom depth is available, in May - Aug
    df = pd.DataFrame(data, index=times[idx])
    df.reset_index(inplace=True)
    df['pH'] = df['pH_TS_insitu_measured']
    df2 = df[df.Depth_bottom != -999]
    df2 = df2[df2.Month_UTC > 4]
    df2 = df2[df2.Month_UTC < 9]
    df2['depth_diff'] = df2.Depth_bottom - df2.Depth
    df2['bottom_threshold'] = df2.Depth_bottom * .8  # keep samples within

    # get rid of data where the sampling depth is deeper than the bottom depth
    df2 = df2[df2['depth_diff'] >= 0]

    dff = df2[df2['depth_diff'] <= 10].copy()

    # use calculated pH if measured isn't available
    for idx, row in dff.iterrows():
        if row['pH'] == -999:
            dff.loc[row.name, 'pH'] = row.pH_TS_insitu_calculated

    dff.replace(-999, np.nan, inplace=True)
    dff.dropna(subset=['pH'], inplace=True)

    # The ECOA dataset doesn't have bottom depth recorded, so have to treat this one differently.
    # get ECOA data
    df_ecoa = df[(df.Cruise_ID == 'ECOA1') | (df.Cruise_ID == 'ECOA2')]

    # use calculated pH if measured isn't available
    for idx, row in df_ecoa.iterrows():
        if row['pH'] == -999:
            df_ecoa.loc[row.name, 'pH'] = row.pH_TS_insitu_calculated

    df_ecoa.replace(-999, np.nan, inplace=True)
    df_ecoa.dropna(subset=['pH'], inplace=True)

    # for each profile, find the deepest sample and compare the depth to the global bathymetry file
    extent = [-78, -65, 35, 45]
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))
    ecoa_data = dict(pressure=np.array([]),
                     lat=np.array([]),
                     lon=np.array([]),
                     pH=np.array([]))
    profile_num = np.unique(df_ecoa.Profile_number)
    for i, pn in enumerate(profile_num):
        df_pn = df_ecoa[df_ecoa.Profile_number == pn]
        maxdepth = np.nanmax(df_pn.Depth)
        if maxdepth > 10:  # sampling depth has to be >10m
            pn_coords = [df_pn[df_pn.Depth == maxdepth].Longitude.values[0],
                         df_pn[df_pn.Depth == maxdepth].Latitude.values[0]]

            # find the closest bathymetry grid location
            lat_idx = abs(bathy.lat.values - pn_coords[1]).argmin()
            lon_idx = abs(bathy.lon.values - pn_coords[0]).argmin()
            bathy_depth = -bathy.elevation[lat_idx, lon_idx].values

            # calculate 20% of the water column from the global bathymetry file
            # if the measured ECOA value is deeper than that threshold, keep the value
            bathy_threshold = bathy_depth - (bathy_depth * .2)
            if maxdepth > bathy_threshold:
                ecoa_data['pressure'] = np.append(ecoa_data['pressure'], maxdepth)
                ecoa_data['lat'] = np.append(ecoa_data['lat'], pn_coords[1])
                ecoa_data['lon'] = np.append(ecoa_data['lon'], pn_coords[0])
                ecoa_data['pH'] = np.append(ecoa_data['pH'], df_pn[df_pn.Depth == maxdepth].pH.values[0])

    # get glider data
    # for each glider, grab the data in the bottom 1m of the profile except when the profile is deeper than a
    # defined threshold (when the glider is off the shelf and no longer sampling the bottom)
    # we're assuming each profile samples the bottom, except off the shelf

    gldata = dict(pressure=np.array([]),
                  lat=np.array([]),
                  lon=np.array([]),
                  pH=np.array([]))
    glider_data = dict(sbu=gldata,
                       ru=gldata)

    for rf in ruf:
        # RU30 - glider rated for 200m, so ignore profiles > 170m (presumably off the shelf)
        ru_data = glider_data['ru']
        ds = xr.open_dataset(rf)
        try:
            profileid = np.unique(ds.profile_id.values)
            profile_identifier = 'profile_id'
            phvar = 'pH'
        except AttributeError:
            profileid = np.unique(ds.profile_time.values)
            profile_identifier = 'profile_time'
            phvar = 'ph_total_shifted'
        for pid in profileid:
            pidx = np.where(ds[profile_identifier].values == pid)[0]
            maxpress = np.nanmax(ds.pressure.values[pidx])
            if np.logical_and(maxpress > 30, maxpress < 170):
                idx = np.where(ds.pressure.values[pidx] > maxpress - 1)[0]
                ru_data['pressure'] = np.append(ru_data['pressure'], maxpress)
                ru_data['lat'] = np.append(ru_data['lat'], np.nanmean(ds.latitude.values[pidx][idx]))
                ru_data['lon'] = np.append(ru_data['lon'], np.nanmean(ds.longitude.values[pidx][idx]))
                ru_data['pH'] = np.append(ru_data['pH'], statistics.median(ds[phvar].values[pidx][idx]))

    for sf in sbuf:
        # SBU - glider rated for 350m but Charlie keeps it at 180m off the shelf, so ignore profiles > 170m
        sbu_data = glider_data['sbu']
        ds = xr.open_dataset(sf)
        ptime = np.unique(ds.profiletime)  # these data are actually grouped by yo (one down-up pair)
        for pt in ptime:
            dst = ds.sel(profiletime=pt)
            maxpress = np.nanmax(dst.pressure.values)
            if np.logical_and(maxpress > 30, maxpress < 170):
                idx = np.where(dst.pressure.values > maxpress - 1)[0]
                sbu_data['pressure'] = np.append(sbu_data['pressure'], maxpress)
                sbu_data['lat'] = np.append(sbu_data['lat'], np.unique(dst.profilelat.values)[0])
                sbu_data['lon'] = np.append(sbu_data['lon'], np.unique(dst.profilelon.values)[0])
                sbu_data['pH'] = np.append(sbu_data['pH'], statistics.median(dst.ph_total.values[idx]))

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

    # plot CODAP data
    sct = ax.scatter(dff.Longitude, dff.Latitude, c=dff.pH, marker='.', s=100, cmap=cmo.cm.matter,
                     transform=ccrs.PlateCarree(), zorder=10)

    # plot ECOA data
    sct = ax.scatter(ecoa_data['lon'], ecoa_data['lat'], c=ecoa_data['pH'], marker='.',
                     s=100, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(), zorder=10)

    # plot glider data
    sct = ax.scatter(glider_data['sbu']['lon'], glider_data['sbu']['lat'], c=glider_data['sbu']['pH'], marker='.',
                     s=75, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(), zorder=10)
    sct = ax.scatter(glider_data['ru']['lon'], glider_data['ru']['lat'], c=glider_data['ru']['pH'], marker='.',
                     s=75, cmap=cmo.cm.matter, transform=ccrs.PlateCarree(), zorder=10)

    # Set colorbar height equal to plot height
    divider = make_axes_locatable(ax)
    cax = divider.new_horizontal(size='5%', pad=0.05, axes_class=plt.Axes)
    fig.add_axes(cax)

    # generate colorbar
    cb = plt.colorbar(sct, cax=cax)
    cb.set_label(label='Bottom pH (summer)')

    sfile = os.path.join(os.path.dirname(f), 'bottom_ph_map-test.png')
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
