#! /usr/bin/env python3

"""
Author: Lori Garzio on 3/8/2021
Last modified: 12/1/2023
"""

import PyCO2SYS as pyco2
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER


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
