#!/usr/bin/env python

"""
Author: Lori Garzio on 10/18/2024
Last modified: 10/6/2026
Filter the CODAP-NA v2026 dataset to select data within the study region, remove questionable and 
bad data, and calculate aragonite saturation state if not available.
CODAP-NA v2026 is available in the NCEI OCADS data portal, Accession 0315529
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
import yaml
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in IDE


def main(lon_bounds, lat_bounds, codap_file):
    # get CODAP data
    codap_file = os.path.join(os.getenv('HOME'), codap_file)
    ds = xr.open_dataset(codap_file)
    
    # create a dataframe of the cruise information
    cruise_mapping = pd.DataFrame(
    {
        "Cruise_key": np.arange(1, len(ds.EXPOCODE_list.values) + 1),
        "EXPOCODE": ds.EXPOCODE_list.values,
        "Cruise_ID": ds.Cruise_ID_list.values
    }
)
    
    ds = ds.drop_dims(['cruise'])
    df = ds.to_dataframe()

    # config file that lists the variables to be included in the filtered dataset
    configfile = os.path.join(os.path.dirname(__file__), 'configs', 'codapv2026vars.yml')
    with open(configfile, "r") as file:
        codap_vars = yaml.safe_load(file)
    codap_vars = codap_vars.split(' ')
    
    # only keep list of variables in the config file and merge the cruise mapping dataframe to get the cruise ID and EXPOCODE
    df = df[codap_vars]
    df2 = df.merge(cruise_mapping, how='outer', on='Cruise_key')

    # filter the dataset to only include data within the study region
    df2['in_study_region'] = False
    for i, row in df2.iterrows():
        coordinate = Point(row.Longitude, row.Latitude)
        if Polygon(list(zip(lon_bounds, lat_bounds))).contains(coordinate):
            df2.loc[row.name, 'in_study_region'] = True

    df = df2.loc[df2["in_study_region"]].drop(columns="in_study_region").copy()
    
    # generate timestamp
    df['year'] = df['Year_UTC'].apply(int)
    df['month'] = df['Month_UTC'].apply(int)
    df['day'] = df['Day_UTC'].apply(int)
    df['time'] = pd.to_datetime(df[['year', 'month', 'day']])

    # flag_values: [2, 3, 6, 9]
    # flag_meanings: good questionable_but_retained questionable missing_or_invalid
    # remove questionable (6) temperature and salinity
    df.loc[df.CTDTEMP_flag == 6, 'CTDTEMP_ITS90_deg_C'] = np.nan
    df.loc[df.recommended_Salinity_flag == 6, 'recommended_Salinity_PSS78'] = np.nan

    # remove questionable (6) pH data then drop rows with no pH data
    # pH_TS_insitu_measured: measured pH adjusted to in situ conditions
    # pH_TS_insitu_calculated: calculated pH adjusted to in situ conditions
    # pH_TS_insitu_combined: pH at in situ conditions from both measured and calculated values
    df.loc[df.pH_TS_insitu_combined_flag == 6, 'pH_TS_insitu_combined'] = np.nan
    df = df.dropna(subset=["pH_TS_insitu_combined"])

    # remove questionable (6) TA data
    df.loc[df.TALK_flag == 6, 'TALK_umol_kg'] = np.nan

    # convert the observation type to a string
    df["Observation_type"] = df["Observation_type"].replace({
        1: "Niskin",
        2: "Flow-through"})

    # If aragonite saturation state isn't available, calculate it
    for idx, row in df.iterrows():
        if pd.isna(row.Aragonite):
            if not pd.isna(row.TALK_umol_kg):
                omega_arag, pco2, revelle = cf.run_co2sys_ta_ph(row.TALK_umol_kg,
                                                                row.pH_TS_insitu_combined,
                                                                row.recommended_Salinity_PSS78,
                                                                row.CTDTEMP_ITS90_deg_C,
                                                                row.CTDPRES_dbar)
                df.loc[row.name, 'Aragonite'] = float(np.round(omega_arag, 2))

    df.to_csv(os.path.join(os.path.dirname(codap_file), 'CODAP-NA_v2026-filtered-for-analysis.csv'), index=False)


if __name__ == '__main__':
    lons = [-78, -65, -65, -78]  # longitude boundaries for grabbing vessel-based data
    lats = [35, 35, 45, 45]  # latitude boundaries for grabbing vessel-based data
    codap = 'rucool/Saba/OA_cruise_data/CODAPv2026/0315529/2.2/data/0-data/SDIS_submission_CODAP_NA_V2026/Output_file.nc'
    main(lons, lats, codap)
