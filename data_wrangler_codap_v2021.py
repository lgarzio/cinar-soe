#!/usr/bin/env python

"""
Author: Lori Garzio on 10/18/2024
Last modified: 6/26/2026
Filter the CODAP-NA v2021 dataset to select data within the study region, remove questionable and 
bad data, and calculate aragonite saturation state if not available.
CODAP-NA v2021 dataset documented here: https://essd.copernicus.org/articles/13/2777/2021/
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in IDE


def main(lon_bounds, lat_bounds, codap_file):
    # get CODAP data
    codap_file = os.path.join(os.getenv('HOME'), codap_file)
    ds = xr.open_dataset(codap_file)
    idx = []
    codap_vars = dict(Day_UTC=np.array([]),
                        Month_UTC=np.array([]),
                        Year_UTC=np.array([]),
                        Cruise_ID=np.array([]),
                        Accession=np.array([]),
                        Observation_type=np.array([]),
                        Profile_number=np.array([]),
                        Latitude=np.array([]),
                        Longitude=np.array([]),
                        CTDPRES=np.array([]),
                        Depth=np.array([]),
                        Depth_bottom=np.array([]),
                        CTDTEMP_ITS90=np.array([]),
                        CTDTEMP_flag=np.array([]),
                        recommended_Salinity_PSS78=np.array([]),
                        recommended_Salinity_flag=np.array([]),
                        pH_TS_insitu_calculated=np.array([]),
                        pH_TS_insitu_measured=np.array([]),
                        pH_flag=np.array([]),
                        fCO2_insitu_calculated_uatm=np.array([]),
                        fCO2_insitu_measured_uatm=np.array([]),
                        fCO2_flag=np.array([]),
                        TALK=np.array([]),
                        TALK_flag=np.array([]),
                        Aragonite=np.array([]))

    # make sure the data are within the defined extent
    for i, lon in enumerate(ds.Longitude.values):
        if Polygon(list(zip(lon_bounds, lat_bounds))).contains(Point(lon, ds.Latitude.values[i])):
            idx.append(i)
            for key in codap_vars.keys():
                if key in ['Cruise_ID', 'Observation_type']:
                    cid = ds[key].values[:, i]
                    cid = [x.decode('UTF-8') for x in cid]
                    codap_vars[key] = np.append(codap_vars[key], ''.join(cid).strip())
                else:
                    codap_vars[key] = np.append(codap_vars[key], ds[key].values[i])

    df = pd.DataFrame(codap_vars)

    # generate timestamp
    df['year'] = df['Year_UTC'].apply(int)
    df['month'] = df['Month_UTC'].apply(int)
    df['day'] = df['Day_UTC'].apply(int)
    df['time'] = pd.to_datetime(df[['year', 'month', 'day']])

    df['pH'] = df['pH_TS_insitu_measured']  # pH recalculated at insitu temperature
    df['pCO2'] = df['fCO2_insitu_measured_uatm']  # pCO2 measured at insitu temperature

    # use calculated pH and pCO2 if measured isn't available
    for idx, row in df.iterrows():
        if row['pH'] == -999:
            df.loc[row.name, 'pH'] = row.pH_TS_insitu_calculated
        if row['pCO2'] == -999:
            df.loc[row.name, 'pCO2'] = row.fCO2_insitu_calculated_uatm

    # remove questionable (3) and bad (4) pH flags
    # remove missing pH values (-999) *can't use the pH flag = 9 because sometimes that's applied when pH values
    # are available
    df = df[df.pH_flag != 3]
    df = df[df.pH_flag != 4]
    df = df[df.pH != -999]

    # remove questionable (3) and bad (4) pCO2 data
    df.loc[df.fCO2_flag == 3, 'pCO2'] = np.nan
    df.loc[df.fCO2_flag == 4, 'pCO2'] = np.nan

    df['aragonite_estimated'] = ''

    # If aragonite saturation state isn't available, calculate it
    for idx, row in df.iterrows():
        if row.Aragonite == -999:
            # WOCE flags: 2 = Acceptable, 3 = Questionable, 6 = Average of duplicates, 9 = Missing
            if np.logical_or(row.TALK_flag == 2, row.TALK_flag == 6):  # if TA flags are acceptable
                if row.pH_flag != 3:  # if pH flags are acceptable
                    omega_arag, pco2, revelle = cf.run_co2sys_ta_ph(row.TALK,
                                                                    row.pH,
                                                                    row.recommended_Salinity_PSS78,
                                                                    row.CTDTEMP_ITS90,
                                                                    row.CTDPRES)
                    df.loc[row.name, 'aragonite_estimated'] = np.round(omega_arag, 2)
                else:
                    df.loc[row.name, 'aragonite_estimated'] = np.nan
            else:
                df.loc[row.name, 'aragonite_estimated'] = np.nan
        else:
            df.loc[row.name, 'aragonite_estimated'] = np.nan

    df.replace(-999, np.nan, inplace=True)
    df.to_csv(os.path.join(os.path.dirname(codap_file), 'CODAP-NA_v2021-filtered-for-analysis.csv'))


if __name__ == '__main__':
    lons = [-78, -65, -65, -78]  # longitude boundaries for grabbing vessel-based data
    lats = [35, 35, 45, 45]  # latitude boundaries for grabbing vessel-based data
    codap = 'rucool/Saba/OA_cruise_data/CODAPv2021/CODAP_NA_v2021.nc'
    main(lons, lats, codap)
