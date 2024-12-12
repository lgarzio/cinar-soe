#! /usr/bin/env python3

"""
Author: Lori Garzio on 3/8/2021
Last modified: 12/12/2024
"""

import PyCO2SYS as pyco2
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.io.shapereader import Reader
from shapely.ops import unary_union


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
