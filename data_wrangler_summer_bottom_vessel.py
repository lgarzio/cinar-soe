#!/usr/bin/env python

"""
Author: Lori Garzio on 10/18/2021
Last modified: 11/17/2023
Grab summer bottom-water pH and omega data from CODAP-NA and additional ECOMON datasets. Export as NetCDF.
CODAP-NA dataset documented here: https://essd.copernicus.org/articles/13/2777/2021/
"""

import os
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
import PyCO2SYS as pyco2
from collections import OrderedDict
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console


def main(lon_bounds, lat_bounds, codap_file, ecomon_files, savedir):
    bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    extent = [-78, -65, 35, 45]
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))

    # initialize dictionary to append bottom pH and aragonite data from cruises (CODAP and ECOMON datasets)
    data = {
        "coords": {
            "time": {"dims": "time", "data": np.array([], dtype='datetime64[ns]')}
        },
        "attrs": {
            "comment": "Synthesis of bottom pH and aragonite saturation state data from vessel-based measurements "
                       "collected during summer months (June-August) that were spatially limited to the U.S. Northeast "
                       "Shelf. ",
            "data_sources": "CODAP-NA dataset (https://essd.copernicus.org/articles/13/2777/2021/) and subsequent "
                            "ECOMON cruises not included in the CODAP-NA dataset. Additional ECOMON data were "
                            "downloaded from NCEI's Ocean Carbon and Acidification Data Portal "
                            "(https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system-portal/)"
        },
        "dims": "time",
        "data_vars": {
            "data_source": {
                "dims": "time",
                "data": np.array([], dtype='<U32'),
                "attrs": {
                    "units": "1",
                    "comment": "Source of data"
                }
            },
            "cruise": {
                "dims": "time",
                "data": np.array([], dtype='<U32'),
                "attrs": {
                    "units": "1",
                    "comment": "Cruise ID from original data source"
                }
            },
            "lat": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "degrees_north",
                    "long_name": "Latitude",
                    "comment": "Latitude from original data source"
                }
            },
            "lon": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "degrees_east",
                    "long_name": "Longitude",
                    "comment": "Longitude from original data source"
                }
            },
            "depth": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "m",
                    "long_name": "Depth",
                    "description": "Depth at which sample was collected",
                    "comment": "From the original data source"
                }
            },
            "bottom_depth": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "m",
                    "long_name": "Bottom Depth",
                    "description": "Depth of water column at the sampling station",
                    "comment": "From the original data source, when provided. When not provided, bottom depth was "
                               "calculated using the sample coordinates and GEBCO’s gridded global bathymetry "
                               "(https://download.gebco.net/)"
                }
            },
            "pH_bottom": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "1",
                    "long_name": "pH",
                    "description": "Bottom pH from original data source",
                    "comment": "Measured in situ pH from the original data source, when available. If measured pH was "
                               "not available, calculated pH was used. For the ECOMON datasets, pH at in situ "
                               "temperature, pressure, and salinity were calculated from pH_TS_20C using PyCO2SYS. "
                               "Bottom pH is the deepest measurement of a "
                               "vertical CTD/Rosette cast where water samples were collected, for profiles deeper "
                               "than 10m, provided the sampling depth was no shallower than the bottom 20% of total "
                               "water column depth."
                }
            },
            "omega_bottom": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "1",
                    "long_name": "Aragonite Saturation State",
                    "description": "Bottom Aragonite Saturation State",
                    "comment": "Bottom aragonite saturation state from the original data source, when available. "
                               "When not available, calculated using PyCO2SYS (Humphreys et al. 2020) with inputs of "
                               "pressure, temperature, salinity, total alkalinity, and pH. "
                               "Bottom aragonite saturation state is the deepest measurement of a "
                               "vertical CTD/Rosette cast where water samples were collected, for profiles deeper "
                               "than 10m, provided the sampling depth was no shallower than the bottom 20% of total "
                               "water column depth."
                }
            },
        },
    }

    # get CODAP data
    ds = xr.open_dataset(codap_file)
    idx = []
    codap_vars = dict(Day_UTC=np.array([]),
                      Month_UTC=np.array([]),
                      Year_UTC=np.array([]),
                      Cruise_ID=np.array([]),
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
                      TALK=np.array([]),
                      TALK_flag=np.array([]),
                      Aragonite=np.array([]))

    # make sure the data are within the defined extent
    for i, lon in enumerate(ds.Longitude.values):
        if Polygon(list(zip(lon_bounds, lat_bounds))).contains(Point(lon, ds.Latitude.values[i])):
            idx.append(i)
            for key in codap_vars.keys():
                if key == 'Cruise_ID':
                    cid = ds.Cruise_ID.values[:, i]
                    cid = [x.decode('UTF-8') for x in cid]
                    codap_vars[key] = np.append(codap_vars[key], ''.join(cid).strip())
                else:
                    codap_vars[key] = np.append(codap_vars[key], ds[key].values[i])

    # select data from June - Aug
    df = pd.DataFrame(codap_vars)
    df = df[df.Month_UTC > 5]
    df = df[df.Month_UTC < 9]
    df['pH'] = df['pH_TS_insitu_measured']

    # generate timestamp
    df['year'] = df['Year_UTC'].apply(int)
    df['month'] = df['Month_UTC'].apply(int)
    df['day'] = df['Day_UTC'].apply(int)
    df['time'] = pd.to_datetime(df[['year', 'month', 'day']])

    # use calculated pH if measured isn't available
    for idx, row in df.iterrows():
        if row['pH'] == -999:
            df.loc[row.name, 'pH'] = row.pH_TS_insitu_calculated

    df = df[df.pH != -999]

    # get rid of questionable pH flags
    df = df[df.pH_flag != 3]

    df['aragonite_final'] = ''

    # If aragonite saturation state is available, use that value. If it's not available, calculate aragonite saturation
    # state
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
                    df.loc[row.name, 'aragonite_final'] = np.round(omega_arag, 2)
                else:
                    df.loc[row.name, 'aragonite_final'] = np.nan
            else:
                df.loc[row.name, 'aragonite_final'] = np.nan
        else:
            df.loc[row.name, 'aragonite_final'] = row.Aragonite

    # for each profile on each cruise, find the deepest sample and compare the sample depth to either the water depth
    # provided in the file, or depth from the global bathymetry file
    cruises = np.unique(df.Cruise_ID)
    for cruise in cruises:
        dfc = df[df.Cruise_ID == cruise]
        profile_num = np.unique(dfc.Profile_number)
        for profile in profile_num:
            dfc_profile = dfc[dfc.Profile_number == profile]
            maxdepth = np.nanmax(dfc_profile.Depth)
            if maxdepth > 10:  # sampling depth has to be >10m
                dfc_profile_max = dfc_profile[dfc_profile.Depth == maxdepth]
                profile_coords = [dfc_profile_max.Longitude.values[0], dfc_profile_max.Latitude.values[0]]
                if dfc_profile_max.Depth_bottom.values[0] != -999:  # compare to the recorded station depth
                    station_water_depth = dfc_profile_max.Depth_bottom.values[0]
                else:  # compare to the global bathymetry file
                    lat_idx = abs(bathy.lat.values - profile_coords[1]).argmin()
                    lon_idx = abs(bathy.lon.values - profile_coords[0]).argmin()
                    station_water_depth = -bathy.elevation[lat_idx, lon_idx].values

                # if the sample is within +/- 20% of the water column, keep the value
                depth_threshold = [station_water_depth * .8, station_water_depth * 1.2]
                if np.logical_and(maxdepth > depth_threshold[0], maxdepth < depth_threshold[1]):
                    tm = dfc_profile_max.time.values[0]
                    pH_bottom = dfc_profile_max.pH.values[0]
                    omega_bottom = np.round(dfc_profile_max.aragonite_final.values[0], 2)
                    data['coords']['time']['data'] = np.append(data['coords']['time']['data'], tm)
                    data['data_vars']['data_source']['data'] = np.append(data['data_vars']['data_source']['data'],
                                                                         'CODAP_NA_v2021')
                    data['data_vars']['cruise']['data'] = np.append(data['data_vars']['cruise']['data'], cruise)
                    data['data_vars']['depth']['data'] = np.append(data['data_vars']['depth']['data'], maxdepth)
                    data['data_vars']['bottom_depth']['data'] = np.append(data['data_vars']['bottom_depth']['data'],
                                                                          station_water_depth)
                    data['data_vars']['lat']['data'] = np.append(data['data_vars']['lat']['data'], profile_coords[1])
                    data['data_vars']['lon']['data'] = np.append(data['data_vars']['lon']['data'], profile_coords[0])
                    data['data_vars']['pH_bottom']['data'] = np.append(data['data_vars']['pH_bottom']['data'],
                                                                       pH_bottom)
                    data['data_vars']['omega_bottom']['data'] = np.append(data['data_vars']['omega_bottom']['data'],
                                                                          omega_bottom)

    # additional EcoMon data that aren't included in CODAP
    for ef in ecomon_files:
        df_ecomon = pd.read_csv(ef)
        df_ecomon = df_ecomon[df_ecomon.Month_UTC > 5]
        df_ecomon = df_ecomon[df_ecomon.Month_UTC < 9]
        df_ecomon.replace(-999, np.nan, inplace=True)
        df_ecomon.dropna(subset=['pH_TS_20C'], inplace=True)

        df_ecomon['year'] = df_ecomon['Year_UTC'].apply(int)
        df_ecomon['month'] = df_ecomon['Month_UTC'].apply(int)
        df_ecomon['day'] = df_ecomon['Day_UTC'].apply(int)
        df_ecomon['time'] = pd.to_datetime(df_ecomon[['year', 'month', 'day']])

        # remove flow-thru data, only use Niskin sampling
        df_ecomon = df_ecomon[df_ecomon.Observation_Type == 'Niskin']

        # combine Station_ID and Cast_number to get a unique profile
        df_ecomon['profile_id'] = df_ecomon.Station_ID.astype(str) + '_' + df_ecomon.Cast_number.astype(str)

        # for each profile, find the deepest sample and compare the sample depth to the water column depth
        # keep the sample if it's within 20% of the bottom depth
        profile_num = np.unique(df_ecomon.profile_id)
        for profile in profile_num:
            df_ecomon_profile = df_ecomon[df_ecomon.profile_id == profile]
            maxdepth = np.nanmax(df_ecomon_profile.Depth_meters)
            if maxdepth > 10:  # sampling depth has to be >10m
                df_ecomon_profile_max = df_ecomon_profile[df_ecomon_profile.Depth_meters == maxdepth]
                profile_coords = [df_ecomon_profile_max.Longitude_Dec_Deg.values[0],
                                  df_ecomon_profile_max.Latitude_Dec_Deg.values[0]]
                station_water_depth = np.nanmax(df_ecomon_profile_max.Depth_Bottom_meters)
                depth_threshold = [station_water_depth * .8, station_water_depth * 1.2]
                if np.logical_and(maxdepth > depth_threshold[0], maxdepth < depth_threshold[1]):
                    tm = df_ecomon_profile_max['time'].values[0]
                    ecomon_cruise = df_ecomon_profile_max.Cruise_ID.values[0]

                    data['coords']['time']['data'] = np.append(data['coords']['time']['data'], tm)
                    data['data_vars']['data_source']['data'] = np.append(data['data_vars']['data_source']['data'],
                                                                         'ECOMON-NCEI')
                    data['data_vars']['cruise']['data'] = np.append(data['data_vars']['cruise']['data'], ecomon_cruise)
                    data['data_vars']['depth']['data'] = np.append(data['data_vars']['depth']['data'], maxdepth)
                    data['data_vars']['bottom_depth']['data'] = np.append(data['data_vars']['bottom_depth']['data'],
                                                                          station_water_depth)
                    data['data_vars']['lat']['data'] = np.append(data['data_vars']['lat']['data'], profile_coords[1])
                    data['data_vars']['lon']['data'] = np.append(data['data_vars']['lon']['data'], profile_coords[0])

                    # calculate pH at in situ temperature, pressure, and salinity using PyCO2SYS
                    ph_20c = df_ecomon_profile_max.pH_TS_20C.values[0]
                    par1_type = 3
                    ta = df_ecomon_profile_max['TA_umol/kg'].values[0]
                    if ta == -999:
                        ta = 2200
                    par2_type = 1

                    kwargs = dict(salinity=df_ecomon_profile_max.CTDSAL_PSS78.values[0],
                                  temperature=20,
                                  temperature_out=df_ecomon_profile_max.CTDTEMP_ITS90.values[0],
                                  pressure=0,
                                  pressure_out=df_ecomon_profile_max.CTDPRES_dbar.values[0],
                                  opt_pH_scale=1,
                                  opt_k_carbonic=4,
                                  opt_k_bisulfate=1,
                                  opt_total_borate=1,
                                  opt_k_fluoride=2)

                    results = pyco2.sys(ph_20c, ta, par1_type, par2_type, **kwargs)
                    pH_bottom = np.round(results['pH_out'], 4)

                    # calculate aragonite
                    omega_bottom, pco2, revelle = cf.run_co2sys_ta_ph(df_ecomon_profile_max['TA_umol/kg'].values[0],
                                                                    pH_bottom,
                                                                    df_ecomon_profile_max.CTDSAL_PSS78.values[0],
                                                                    df_ecomon_profile_max.CTDTEMP_ITS90.values[0],
                                                                    df_ecomon_profile_max.CTDPRES_dbar.values[0])

                    data['data_vars']['pH_bottom']['data'] = np.append(data['data_vars']['pH_bottom']['data'],
                                                                       pH_bottom)
                    data['data_vars']['omega_bottom']['data'] = np.append(data['data_vars']['omega_bottom']['data'],
                                                                          np.round(omega_bottom, 2))

    # save as netcdf
    outds = xr.Dataset.from_dict(data)

    # add created time to global attrs
    datetime_format = '%Y-%m-%dT%H:%M:%SZ'
    created = dt.datetime.utcnow().strftime(datetime_format)  # creation time Timestamp
    time_start = pd.to_datetime(np.nanmin(outds.time.values)).strftime(datetime_format)
    time_end = pd.to_datetime(np.nanmax(outds.time.values)).strftime(datetime_format)

    global_attributes = OrderedDict([
        ('date_created', created),
        ('date_modified', created),
        ('time_coverage_start', time_start),
        ('time_coverage_end', time_end),
        ('creator_email', 'lgarzio@marine.rutgers.edu'),
        ('creator_name', 'Lori Garzio'),
        ('creator_url', 'rucool.marine.rutgers.edu'),
        ('institution', 'Rutgers University'),
        ('contributor_name', 'Grace Saba,Lori Garzio'),
        ('contributor_role', 'Principal Investigator,Data Management')
    ])

    global_attributes.update(outds.attrs)

    outds = outds.assign_attrs(global_attributes)
    outds = outds.sortby(outds.time)

    # Add compression to all variables
    encoding = {}
    for k in outds.data_vars:
        encoding[k] = {'zlib': True, 'complevel': 1}

    encoding['time'] = dict(units='seconds since 1970-01-01 00:00:00', calendar='gregorian', zlib=False,
                            _FillValue=False, dtype=np.double)

    save_file = os.path.join(savedir, 'vessel_based_summer_bottom_OA_data.nc')
    outds.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims='time')


if __name__ == '__main__':
    lons = [-78, -65, -65, -78]  # longitude boundaries for grabbing vessel-based data
    lats = [35, 35, 45, 45]  # latitude boundaries for grabbing vessel-based data
    codap = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/CODAP_NA_v2021.nc'
    ecomon = ['/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/EcoMon/33GG20190815-GU1902_data.csv',
              '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/EcoMon/33HH20190522-HB1902_data.csv',
              '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/EcoMon/334B20210805-PC2104_Data.csv',
              '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/EcoMon/33HH20220531_HB2204_Data.csv']
    save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/output_nc'
    main(lons, lats, codap, ecomon, save_directory)
