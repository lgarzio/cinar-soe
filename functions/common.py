#! /usr/bin/env python3

"""
Author: Lori Garzio on 3/8/2021
Last modified: 6/26/2026
"""

import PyCO2SYS as pyco2
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.io.shapereader import Reader
from shapely.ops import unary_union
import numpy as np


def add_map_features(axis, extent, edgecolor=None, oceancolor='none'):
    edgecolor = edgecolor or 'black'

    land = cfeature.NaturalEarthFeature('physical', 'land', '10m', edgecolor=edgecolor, facecolor='tan')

    state_lines = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='10m',
        facecolor='none')

    # Axes properties and features
    axis.set_extent(extent)
    axis.set_facecolor(oceancolor)
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


def calc_ta_nyb(season, sal):
    # calculate TA for the New York Bight region
    if season == 'DJF':
        m = 48.71
        b = 601.75
    elif season == 'MAM':
        m = 42.37
        b = 817.59
    elif season == 'JJA':
        m = 50.75
        b = 541.81
    elif season == 'SON':
        m = 46.42
        b = 687.99

    ta = m * sal + b
    return ta


def calc_ta_gom(temp, sal):
    # calculate TA for the Gulf of Maine region according to McGarry et al 2020
    # normalize temperature and salinity to the McGarry et al 2020 data from Table 3
    tempn = (temp - 13.20) / 5.92
    saln = (sal - 34.40) / 1.49
    ta = 2289 + (0.758 * tempn) + (69.2 * saln)

    return ta


def return_noaa_polygons():
    """
    Combine multiple strata levels defined in the NOAA bottom trawl survey strata downloaded from
    https://github.com/NOAA-EDAB/FisheryConditionLinks/tree/master/NES_BOTTOM_TRAWL_STRATA
    for the NYB region into inshore, midshelf, and offshore. STRATA mapping from Laura Nazzaro.
    """
    shpfile = Reader(
        '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/NES_BOTTOM_TRAWL_STRATA/NES_BOTTOM_TRAWL_STRATA.shp')
    shplist = list(shpfile.records())

    strata_mapping = dict(
        inshore=dict(
            snum=[3150, 3160, 3170, 3180, 3190, 3200, 3120, 3130, 3140, 3100, 3090, 3110, 3060, 3070, 3080],
            color='#9e36d7ff',
            poly=[]
        ),  # purple
        midshelf=dict(snum=[1730, 1010], color='#d7369eff', poly=[]),  # pink
        offshore=dict(snum=[1740, 1750, 1760, 1020, 1030, 1040], color='#95c983ff', poly=[])  # green
    )

    # combine regions to NYB inshore, midshelf, offshore
    for region_code, sm in strata_mapping.items():
        polys = []
        for sl in shplist:
            if sl.attributes['STRATA'] in sm['snum']:
                poly = sl.geometry
                polys.append(poly)

        outside_poly = unary_union(polys)
        strata_mapping[region_code]['poly'] = outside_poly

    return strata_mapping


def return_bottom_data(df,
                       bathymetry_file,
                       profile_coords,
                       depth_col='Depth',
                       bottom_depth_col='Depth_bottom',
                       pH_col='pH', 
                       omega_col='Aragonite', 
                       omega_est_col='aragonite_estimated', 
                       temp_col='CTDTEMP_ITS90'):
    """
    Returns data points collected at the bottom of a profile
    :param df: dataframe containing data for one profile
    :param bathymetry_file: xarray dataset containing global bathymetry data
    :param profile_coords: tuple containing the coordinates of the profile (lon, lat)
    :param depth_col: name of the column containing depth data, default='Depth'
    :param bottom_depth_col: name of the column containing bottom depth data, default='Depth_bottom'
    :param pH_col: name of the column containing pH data, default='pH'
    :param omega_col: name of the column containing omega data, default='Aragonite'
    :param omega_est_col: name of the column containing estimated omega data, default='aragonite_estimated'
    :param temp_col: name of the column containing temperature data, default='CTDTEMP_ITS90'
    """
    max_depth = np.nanmax(df[depth_col])
    if max_depth > 10:  # sampling depth has to be >10m
        dfc = df[df[depth_col] == max_depth]
        try:
            if dfc[bottom_depth_col].values[0] != -999:  # compare to the recorded station depth
                station_water_depth = dfc[bottom_depth_col].values[0]
            else:  # compare to the global bathymetry file
                lat_idx = abs(bathymetry_file.lat.values - profile_coords[1]).argmin()
                lon_idx = abs(bathymetry_file.lon.values - profile_coords[0]).argmin()
                station_water_depth = -bathymetry_file.elevation[lat_idx, lon_idx].values
        except KeyError:
            lat_idx = abs(bathymetry_file.lat.values - profile_coords[1]).argmin()
            lon_idx = abs(bathymetry_file.lon.values - profile_coords[0]).argmin()
            station_water_depth = -bathymetry_file.elevation[lat_idx, lon_idx].values

        # if the sample is within +/- 20% of the water column, keep the value
        depth_threshold = [station_water_depth * .8, station_water_depth * 1.2]
        if np.logical_and(max_depth > depth_threshold[0], max_depth < depth_threshold[1]):
            if len(dfc) > 1:
                # drop lines where measured omega isn't available
                dfc = dfc[dfc[omega_col] > 0]
                if len(dfc) < 1:
                    raise(ValueError)

            pH_bottom = np.nanmedian(np.array(dfc[pH_col]))

            omega_bottom = np.nanmedian(np.array(dfc[omega_col]))

            # if measured omega isn't available (-999) use estimated aragonite
            if bool(omega_bottom < 0):
                omega_bottom = np.nanmedian(np.array(dfc[omega_est_col]))

            temp_bottom = np.nanmedian(np.array(dfc[temp_col]))
            if bool(temp_bottom < 0):
                temp_bottom = np.nan

        else:
            pH_bottom = np.nan
            omega_bottom = np.nan
            temp_bottom = np.nan

    else:
        station_water_depth = np.nan
        pH_bottom = np.nan
        omega_bottom = np.nan
        temp_bottom = np.nan

    return max_depth, pH_bottom, omega_bottom, temp_bottom, station_water_depth


def return_surface_data(df, 
                        depth_col='Depth', 
                        pH_col='pH', 
                        omega_col='Aragonite', 
                        omega_est_col='aragonite_estimated',
                        temp_col='CTDTEMP_ITS90'):
    """
    Returns data points collected at the surface of a profile
    :param df: dataframe containing data for one profile
    :param depth_col: name of the column containing depth data, default='Depth'
    :param pH_col: name of the column containing pH data, default='pH'
    :param omega_col: name of the column containing omega data, default='Aragonite'
    :param omega_est_col: name of the column containing estimated omega data, default='aragonite_estimated'
    """
    min_depth = np.nanmin(df[depth_col])
    if bool(min_depth < 10):
        dfc = df[df[depth_col]== min_depth]

        # drop lines where measured omega isn't available
        if len(dfc) > 1:
            print('check')
            dfc = dfc[dfc[omega_col] > 0]

            # if you removed all rows of data, go back to the original
            if len(dfc) < 1:
                dfc = df[df[depth_col] == min_depth]

        pH_surface = np.nanmedian(np.array(dfc[pH_col]))

        omega_surface = np.nanmedian(np.array(dfc[omega_col]))

        # if measured omega isn't available (-999) use estimated aragonite
        if bool(omega_surface < 0):
            omega_surface = np.nanmedian(np.array(dfc[omega_est_col]))

        temp_surface = np.nanmedian(np.array(dfc[temp_col]))
        if bool(temp_surface < 0):
            temp_surface = np.nan
    else:
        pH_surface = np.nan
        omega_surface = np.nan
        temp_surface = np.nan
    
    return min_depth, pH_surface, omega_surface, temp_surface


def run_co2sys_ta_ph(ta, ph, sal, temp=25, press_dbar=0):
    """
    Runs the PyCO2SYS function using input TA and pH data.
    opt_pH_scale=1 is the default (Total scale), including here for clarity.
    opt_k_carbonic=4 is MCHP73 (Mehrbach et al 1973) refit by DM87 (Dickson & Millero 1987)
    opt_k_bisulfate=1 is the default (Dickson 1990), including here for clarity
    opt_total_borate=1 is the default (Uppstrom 1974), including here for clarity
    opt_k_fluoride=2 is PF87 (Perez & Fraga 1987)
    :param ta: Total Alkalinity array
    :param ph: pH array
    :param sal: salinity array
    :param temp: temperature array, default=25
    :param press_dbar: pressure (in units=dbar) array, default=0
    :return: calculated aragonite saturation, pCO2, revelle factor
    """
    # define input conditions
    par1 = ta  # Total Alkalinity
    par1_type = 1  # parameter 1 type (TA)
    par2 = ph
    par2_type = 3  # parameter 2 type (pH)

    kwargs = dict(salinity=sal,
                  temperature=temp,
                  pressure=press_dbar,
                  opt_pH_scale=1,
                  opt_k_carbonic=4,
                  opt_k_bisulfate=1,
                  opt_total_borate=1,
                  opt_k_fluoride=2)

    results = pyco2.sys(par1, par2, par1_type, par2_type, **kwargs)
    omega_arag = results['saturation_aragonite']  # aragonite saturation state
    pco2 = results['pCO2']  # units = uatm
    revelle = results['revelle_factor']

    return omega_arag, pco2, revelle
