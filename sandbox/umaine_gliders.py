#!/usr/bin/env python

import os
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
import pytz
import json
import netCDF4 as nc
pd.set_option('display.width', 320, "display.max_columns", 10)  # for display in pycharm console


def phcalc(vrs, press, temp, salt, k0, k2, pcoefs):
    """
    Translated from Seabird code for calculating pH from glider outputs
    Translated by: Lori Garzio on 3/3/2021
    Last modified: 3/3/2021
    :param vrs: voltage between reference electrode and ISFET source
    :param press: pressure in decibars
    :param temp: temperature in degrees C
    :param salt: salinity (usually CTD salinity on the PSS)
    :param k0: sensor reference potential (intercept at Temp = 0 C)
    :param k2: linear temperature coefficient (slope)
    :param pcoefs: sensor dependent pressure coefficients
    :return: phfree, phtot
    """
    # SET CONSTANTS
    # Universal gas constant, (R), http://physics.nist.gov/cgi-bin/cuu/Value?r
    R = 8.31446  # J/(mol K)
    F = 96485  # Faraday constant Coulomb / mol
    Tk = 273.15 + temp  # convert to degrees Kelvin
    ln10 = np.log(10)  # natural log of 10

    # CALCULATE PHYSICAL AND THERMODYNAMIC DATA
    # Dickson, A. G., Sabine, C. L., & Christian, J. R. (2007). Guide to best practices for ocean CO2 measurements.

    # IONIC STRENGTH OF SEAWATER (mol / kg H2O)
    # Verified units by comparing to Dickson et al. 2007: Chap 5, p10 Table 2
    # Dickson et al. 2007: Chap 5, p13 Eq 34
    IonS = 19.924 * salt / (1000 - 1.005 * salt)

    # MEAN SEAWATER SULFATE CONCENTRATION (mol / kg solution)
    # This wants to be mol/kg sw  as KHSO4 is on that scale
    # Dickson et al. 2007: Chap 5, p10 Table 2
    Stotal = (0.14 / 96.062) * (salt / 1.80655)

    # MEAN SEAWATER CHLORIDE CONCENTRATION  (mol / kg H20)
    # this wants to be mol/kg H2O as activity is on mol/kg H2O scale
    # Dickson et al. 2007: Chap 5, p10 Table 2
    Cltotal = 0.99889 / 35.453 * salt / 1.80655  # (mol / kg solution)
    Cltotal = Cltotal / (1 - 0.001005 * salt)  # (mol / kg H20)

    # BISULFIDE DISSCIATION CONSTANT AT T,S AND IONIC STRENGTH(mol/kg solution)
    # Dickson et al. 2007: Chap 5, p12 Eq 33
    Khso4 = np.exp(-4276.1 / Tk + 141.328 - 23.093 * np.log(Tk) + (-13856. / Tk + 324.57 - 47.986 * np.log(Tk))
                   * IonS**0.5 + (35474. / Tk - 771.54 + 114.723 * np.log(Tk)) * IonS - 2698 / Tk * IonS**1.5
                   + 1776 / Tk * IonS**2 + np.log(1 - 0.001005 * salt))

    # Millero 1983 Chemical Oceanography vol 8
    # partial molar volume and compressibility of HSO4 in seawater.   delta v is cm^3
    deltaVHSO4 = -18.03 + 0.0466 * temp + 0.000316 * temp**2
    KappaHSO4 = (-4.53 + 0.09 * temp) / 1000
    lnKhso4fac = (-deltaVHSO4 + 0.5 * KappaHSO4 * (press / 10)) * (press / 10) / (R * 10. * Tk)

    # bisulfate association constant at T, S, P
    Khso4TPS = Khso4 * np.exp(lnKhso4fac)

    # GAMMA +/- HCl, activity coefficient of HCl at T/S, P=1
    # ADH is the Debye Huckel constant, calcualted as a polynomial fit to data in Khoo et al. 1977, doi:10.1021/ac50009a016
    # See Martz et al. 2010, DOI 10.4319/lom.2010.8.172, p175
    # Typo in paper 2nd term should be e-4 not e-6
    ADH = (3.4286e-6 * temp**2 + 6.7524e-4 * temp + 0.49172143)
    log10gammaHCl = -ADH * np.sqrt(IonS) / (1 + 1.394 * np.sqrt(IonS)) + (0.08885 - 0.000111 * temp) * IonS

    # Millero 1983 partial molar volume of HCl in seawater. deltaV is cm^3
    deltaVHCl = 17.85 + 0.1044 * temp - 0.001316 * temp**2

    # make activity coefficient of HCl more accurate by including the effect of pressure on partial molar volume of HCl
    # last divide by 10 is for units in cm^3 vs m^3 and the Pascal vs bar units in mks constants. (Meter/kilo/second)
    log10gammaHCLtP = log10gammaHCl + deltaVHCl * (press / 10) / (R * Tk * ln10) / 2. / 10

    # Sensor reference potential
    k0T = k0 + k2 * temp  # temp in deg C

    # CALCULATE PRESSURE CORRECTION (POLYNOMIAL FUNCTION OF PRESSURE)
    # ALL SENSORS HAVE A PRESSURE RESPONSE WHICH IS DETERMINED IN THE LAB
    # AND CONTAINED IN THE POLYNOMIAL Pcoefs
    # pc = [flipud(Pcoefs);0]; % Matlab wants descending powers & n+1 (add 0)
    # pcorr = polyval(pc,Press)
    k0TP = k0T + pcoefs

    # pH on free scale then corrected to get to pH total on mol/kg sw scale
    # pHinsituFree = (Vrs - k0TP) / (R * Tk / F * ln10) + log(Cltotal) / ln10 + 2 * log10gammaHCLtP
    # this will be mol kg H2O. need to convert to mol/kg sw
    phfree = (vrs - k0TP) / (R * Tk / F * ln10) + np.log(Cltotal) / ln10 + 2 * log10gammaHCLtP  # mol/kg-H2O scale

    # CONVERT TO mol/kg-sw scale - JP 2/4/16
    phfree = phfree - np.log10(1 - 0.001005 * salt)  # mol/kg-sw scale

    # convert to total proton scale
    phtot = phfree - np.log10(1 + Stotal / Khso4TPS)

    return phfree, phtot


f = '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/um_242-20210630T1916/delayed/um_242-20210630T1916-profile-sci-delayed.nc'
calfile = '/Users/garzio/Documents/repo/lgarzio/phglider/calibration/sbe10490_20200818.txt'
with open(calfile) as json_file:
    cc = json.load(json_file)

ds = xr.open_dataset(f)
ds = ds.drop_vars(names=['profile_id', 'rowSize'])
ds = ds.swap_dims({'obs': 'time'})
ds = ds.sortby(ds.time)

# convert pH voltages of 0.0 to nan
ds['sbe41n_ph_ref_voltage'][ds['sbe41n_ph_ref_voltage'] == 0.0] = np.nan

df = ds.to_dataframe()

# convert profile_time to seconds since 1970-01-01 so you don't lose the column in resampling
df['profile_time'] = df['profile_time'].astype('int64')//1e9
df = df.resample('1s').mean()

col_subset = ['conductivity', 'salinity', 'sci_water_pressure', 'temperature', 'sbe41n_ph_ref_voltage']
# drop timestamps with duplicated and missing data
df = df.dropna(axis=0, how='all', subset=col_subset)
df.drop_duplicates(inplace=True)

# calculate pH and add to dataframe (pressure is dbar)
df['f_p'] = np.polyval([cc['f6'], cc['f5'], cc['f4'], cc['f3'], cc['f2'], cc['f1'], 0], df.pressure)
phfree, phtot = phcalc(df.sbe41n_ph_ref_voltage, df.pressure, df.temperature, df.salinity, cc['k0'], cc['k2'], df.f_p)
df['ph_total'] = phtot

# drop rows with invalid pH values
df = df[(df.ph_total > 6.5) & (df.ph_total < 9)]

# convert the dataframe to an xarray dataset
ds_mod = df.to_xarray()
ds_mod.to_netcdf('{}_calculated.nc'.format(f.split('.')[0]))