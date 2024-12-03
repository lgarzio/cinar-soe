#!/usr/bin/env python

"""
Author: Lori Garzio on 10/18/2024
Last modified: 10/23/2024
Grab vessel-based bottom-water pH and omega data from CODAP-NA and additional ECOMON and ECOA datasets. Export as NetCDF.
CODAP-NA v2021 dataset documented here: https://essd.copernicus.org/articles/13/2777/2021/
Additional cruise datasets were downloaded from the NCEI OCADs data portal
(https://www.ncei.noaa.gov/products/ocean-carbon-acidification-data-system)
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
from gsw import p_from_z
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console


def main(lon_bounds, lat_bounds, codap_file, extra_files, savedir):
    bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    extent = [-78, -65, 35, 45]
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))

    # initialize dictionary to append bottom pH and aragonite data from cruises (CODAPv2021, ECOMON and ECOA datasets)
    data = {
        "coords": {
            "time": {"dims": "time",
                     "data": np.array([], dtype='float32'),
                     "attrs": {
                         "units": "seconds since 1970-01-01T00:00:00Z",
                         "time_origin": "01-JAN-1970 00:00:00"
                        }
                     }
        },
        "attrs": {
            "comment": "Synthesis of bottom pH and aragonite saturation state data from vessel-based measurements "
                       "that were spatially limited to the U.S. Northeast Shelf.",
            "data_sources": "CODAP-NA v2021 dataset (https://essd.copernicus.org/articles/13/2777/2021/) and subsequent "
                            "ECOMON and ECOA cruises not included in the CODAP-NA v2021 dataset. Additional datasets "
                            "were downloaded from NCEI's Ocean Carbon and Acidification Data Portal "
                            "(https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system-portal/)"
        },
        "dims": "time",
        "data_vars": {
            "data_source": {
                "dims": "time",
                "data": np.array([], dtype='object'),
                "attrs": {
                    "units": "1",
                    "comment": "Source of data"
                }
            },
            "cruise": {
                "dims": "time",
                "data": np.array([], dtype='object'),
                "attrs": {
                    "units": "1",
                    "comment": "Cruise ID"
                }
            },
            "obs_type": {
                "dims": "time",
                "data": np.array([], dtype='object'),
                "attrs": {
                    "units": "1",
                    "comment": "Observation type"
                }
            },
            "accession": {
                "dims": "time",
                "data": np.array([], dtype='int32'),
                "attrs": {
                    "units": "1",
                    "comment": "NCEI Accession number"
                }
            },
            "lat": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "degrees_north",
                    "long_name": "Latitude"
                }
            },
            "lon": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "degrees_east",
                    "long_name": "Longitude"
                }
            },
            "depth": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "m",
                    "long_name": "Depth",
                    "description": "Depth at which sample was collected"
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
                               "determined using the sample coordinates and GEBCO’s gridded global bathymetry "
                               "(https://download.gebco.net/)"
                }
            },
            "pH_bottom": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "1",
                    "long_name": "Bottom pH",
                    "comment": "For the CODAP-NA v2021 dataset, pH_TS_insitu_measured was used when available. "
                               "Otherwise, pH_TS_insitu_calculated was used. For the datasets downloaded from "
                               "NCEI, pH at in situ temperature, pressure, and salinity were calculated from "
                               "pH_TS_20C using PyCO2SYS.",
                    "description": "Bottom pH is the deepest measurement of a vertical CTD/Rosette cast where water "
                                   "samples were collected for profiles deeper than 10m, provided the sampling depth "
                                   "was no shallower than the bottom 20% of total water column depth. For stations "
                                   "with multiple samples at the deepest measurement, the median value was calculated.",
                }
            },
            "omega_bottom": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "1",
                    "long_name": "Bottom Aragonite Saturation State",
                    "description": "Bottom aragonite saturation state is the deepest measurement of a vertical "
                                   "CTD/Rosette cast where water samples were collected for profiles deeper than 10m, "
                                   "provided the sampling depth was no shallower than the bottom 20% of total "
                                   "water column depth. For stations with multiple samples at the deepest measurement, "
                                   "the median value was calculated.",
                    "comment": "Bottom aragonite saturation state from the original data source, when available. "
                               "When not available, calculated using PyCO2SYS (Humphreys et al. 2020) with inputs of "
                               "pressure, temperature, salinity, total alkalinity, and pH. "

                }
            },
            "temperature_bottom": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "C",
                    "long_name": "Bottom Temperature",
                    "description": "Bottom seawater temperature is the deepest measurement of a vertical "
                                   "CTD/Rosette cast where water samples were collected for profiles deeper than 10m, "
                                   "provided the sampling depth was no shallower than the bottom 20% of total "
                                   "water column depth. For stations with multiple samples at the deepest measurement, "
                                   "the median value was calculated."
                }
            },
        },
    }

#######################################################################################################################
    # get CODAP data
    try:
        # see if the summary csv file already exists
        df = pd.read_csv(os.path.join(os.path.dirname(codap_file), 'CODAP-NA_v2021-filtered-for-bottom-analysis.csv'))
    except FileNotFoundError:
        # if not, create it
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

        # remove flow-thru data, only use Niskin sampling
        df = df[df.Observation_type == 'Niskin']

        df['pH'] = df['pH_TS_insitu_measured']  # pH recalculated at insitu temperature

        # generate timestamp
        df['year'] = df['Year_UTC'].apply(int)
        df['month'] = df['Month_UTC'].apply(int)
        df['day'] = df['Day_UTC'].apply(int)
        df['time'] = pd.to_datetime(df[['year', 'month', 'day']])

        # use calculated pH if measured isn't available
        for idx, row in df.iterrows():
            if row['pH'] == -999:
                df.loc[row.name, 'pH'] = row.pH_TS_insitu_calculated

        # remove questionable (3) and bad (4) pH flags
        # remove missing pH values (-999) *can't use the pH flag = 9 because sometimes that's applied when pH values
        # are available
        df = df[df.pH_flag != 3]
        df = df[df.pH_flag != 4]
        df = df[df.pH != -999]

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

        df.to_csv(os.path.join(os.path.dirname(codap_file), 'CODAP-NA_v2021-filtered-for-bottom-analysis.csv'))

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
                    if len(dfc_profile_max) > 1:
                        # drop lines where measured omega isn't available
                        dfc_profile_max = dfc_profile_max[dfc_profile_max.Aragonite > 0]
                        if len(dfc_profile_max) < 1:
                            raise(ValueError)
                    tm = pd.to_datetime(dfc_profile_max.time.values[0]).timestamp()
                    pH_bottom = np.nanmedian(np.array(dfc_profile_max.pH))

                    omega_bottom = np.nanmedian(np.array(dfc_profile_max.Aragonite))

                    # if measured omega isn't available (-999) use estimated aragonite
                    if bool(omega_bottom < 0):
                        omega_bottom = np.nanmedian(np.array(dfc_profile_max.aragonite_estimated))

                    temp_bottom = np.nanmedian(np.array(dfc_profile_max.CTDTEMP_ITS90))
                    if bool(temp_bottom < 0):
                        temp_bottom = np.nan

                    # add data to dictionary
                    data['coords']['time']['data'] = np.append(data['coords']['time']['data'], tm)

                    data['data_vars']['data_source']['data'] = np.append(data['data_vars']['data_source']['data'],
                                                                         'CODAP_NA_v2021')
                    data['data_vars']['cruise']['data'] = np.append(data['data_vars']['cruise']['data'], cruise)
                    data['data_vars']['obs_type']['data'] = np.append(data['data_vars']['obs_type']['data'],
                                                                      dfc_profile_max.Observation_type.values[0])
                    # data['data_vars']['data_source']['data'].append('CODAP_NA_v2021')
                    # data['data_vars']['cruise']['data'].append(cruise)
                    # data['data_vars']['obs_type']['data'].append(dfc_profile_max.Observation_type.values[0])
                    data['data_vars']['accession']['data'] = np.append(data['data_vars']['accession']['data'],
                                                                       int(dfc_profile_max.Accession.values[0]))
                    data['data_vars']['depth']['data'] = np.append(data['data_vars']['depth']['data'], maxdepth)
                    data['data_vars']['bottom_depth']['data'] = np.append(data['data_vars']['bottom_depth']['data'],
                                                                          station_water_depth)
                    data['data_vars']['lat']['data'] = np.append(data['data_vars']['lat']['data'], profile_coords[1])
                    data['data_vars']['lon']['data'] = np.append(data['data_vars']['lon']['data'], profile_coords[0])
                    data['data_vars']['pH_bottom']['data'] = np.append(data['data_vars']['pH_bottom']['data'],
                                                                       pH_bottom)
                    data['data_vars']['omega_bottom']['data'] = np.append(data['data_vars']['omega_bottom']['data'],
                                                                          omega_bottom)
                    data['data_vars']['temperature_bottom']['data'] = np.append(
                        data['data_vars']['temperature_bottom']['data'], temp_bottom)

#######################################################################################################################
    # add extra datasets
    accession_mapping = {'HB1902': 209045,
                         'GU1902': 209156,
                         'GU1905': 210238,
                         'GU2102': 248269,
                         'PC2104': 249432,
                         'PC2106': 249517,
                         'PC2205': 283758,
                         'HB2204': 276023,
                         'ECOA3': 283329,
                         'HB2302': 296717}

    # additional datasets that aren't included in CODAP v2021
    for ef in extra_files:
        print(ef)
        df = pd.read_csv(ef)
        df.replace(-999, np.nan, inplace=True)
        try:
            df.dropna(subset=['pH_TS_20C'], inplace=True)
            phvar = 'pH_TS_20C'
            tavar = 'TA_umol/kg'
            report_temp = 20
        except KeyError:
            try:
                df.dropna(subset=['pH 20C'], inplace=True)
                phvar = 'pH 20C'
                tavar = 'TA (umol/kg)'
                report_temp = 20
            except KeyError:
                # ECOA-3 dataset
                # ECOA2 and ECOA3 - pH measurements were made at 25C (pers comm Shawn Shellito <Shawn.Shellito@unh.edu>)
                df.dropna(subset=['pH_T_measured_25C'], inplace=True)
                # get rid of questionable (3) and bad (4) pH flags
                df = df[df.pH_flag != 3]
                df = df[df.pH_flag != 4]

                phvar = 'pH_T_measured_25C'
                tavar = 'TA'
                report_temp = 25

                # calculate pressure from depth
                df['CTDPRES_dbar'] = p_from_z(-df['Depth'], df['Latitude'])

        # format date
        try:
            df['year'] = df['Year_UTC'].apply(int)
            df['month'] = df['Month_UTC'].apply(int)
            df['day'] = df['Day_UTC'].apply(int)
            df['time'] = pd.to_datetime(df[['year', 'month', 'day']])
        except KeyError:
            df['time'] = pd.to_datetime(df['Date_UTC'], format="%m/%d/%y")

        # remove flow-thru data, only use Niskin sampling
        df = df[df.Observation_Type == 'Niskin']

        # make sure the data are within the defined extent
        df['in_region'] = ''
        for i, row in df.iterrows():
            try:
                lon = row.Longitude_Dec_Deg
                lat = row.Latitude_Dec_Deg
            except AttributeError:
                lon = row.Longitude
                lat = row.Latitude
            if Polygon(list(zip(lon_bounds, lat_bounds))).contains(Point(lon, lat)):
                df.loc[i, 'in_region'] = 'yes'
            else:
                df.loc[i, 'in_region'] = 'no'

        # drop data if it's not in the region specified
        df = df[df.in_region == 'yes']

        # combine Station_ID and Cast_number to get a unique profile
        try:
            df['profile_id'] = df.Station_ID.astype(str) + '_' + df.Cast_number.astype(str)
        except AttributeError:
            df['profile_id'] = df.STNNBR.astype(str) + '_' + df.CASTNO.astype(str)

        # for each profile, find the deepest sample and compare the sample depth to the water column depth
        # keep the sample if it's within 20% of the bottom depth
        profile_num = np.unique(df.profile_id)
        for profile in profile_num:
            df_profile = df[df.profile_id == profile]

            # find the column name for sample depth
            sample_depth_vars = ['Depth_meters', 'Depth_sampling (M)', 'Depth']
            for sample_depth_var in sample_depth_vars:
                try:
                    maxdepth = np.nanmax(df_profile[sample_depth_var])
                    depthvar = sample_depth_var
                except KeyError:
                    continue

            if maxdepth > 10:  # sampling depth has to be >10m
                df_profile_max = df_profile[df_profile[depthvar] == maxdepth]

                try:
                    profile_coords = [df_profile_max.Longitude_Dec_Deg.values[0],
                                      df_profile_max.Latitude_Dec_Deg.values[0]]
                    station_water_depth = np.nanmax(df_profile_max.Depth_Bottom_meters)
                except AttributeError:
                    try:
                        profile_coords = [df_profile_max.Longitude.values[0],
                                          df_profile_max.Latitude.values[0]]
                        station_water_depth = np.nanmax(df_profile_max['Depth_station (M)'])
                    except KeyError:
                        # ECOA-3 dataset doesn't include station water depth, so grab it from the global bathymetry file
                        profile_coords = [df_profile_max.Longitude.values[0],
                                          df_profile_max.Latitude.values[0]]
                        lat_idx = abs(bathy.lat.values - profile_coords[1]).argmin()
                        lon_idx = abs(bathy.lon.values - profile_coords[0]).argmin()
                        station_water_depth = -bathy.elevation[lat_idx, lon_idx].values

                depth_threshold = [station_water_depth * .8, station_water_depth * 1.2]
                if np.logical_and(maxdepth > depth_threshold[0], maxdepth < depth_threshold[1]):
                    tm = pd.to_datetime(df_profile_max['time'].values[0]).timestamp()
                    cruise = df_profile_max.Cruise_ID.values[0]

                    data['coords']['time']['data'] = np.append(data['coords']['time']['data'], tm)
                    data['data_vars']['data_source']['data'] = np.append(data['data_vars']['data_source']['data'],
                                                                         'NCEI-OCADS')
                    data['data_vars']['cruise']['data'] = np.append(data['data_vars']['cruise']['data'], cruise)
                    data['data_vars']['obs_type']['data'] = np.append(data['data_vars']['obs_type']['data'],
                                                                      df_profile_max.Observation_Type.values[0])
                    cruise_accession = accession_mapping[cruise]
                    data['data_vars']['accession']['data'] = np.append(data['data_vars']['accession']['data'], cruise_accession)

                    data['data_vars']['depth']['data'] = np.append(data['data_vars']['depth']['data'], maxdepth)
                    data['data_vars']['bottom_depth']['data'] = np.append(data['data_vars']['bottom_depth']['data'],
                                                                          station_water_depth)
                    data['data_vars']['lat']['data'] = np.append(data['data_vars']['lat']['data'], profile_coords[1])
                    data['data_vars']['lon']['data'] = np.append(data['data_vars']['lon']['data'], profile_coords[0])

                    # calculate pH at in situ temperature, pressure, and salinity using PyCO2SYS
                    ph_20c = df_profile_max[phvar].values[0]
                    par1_type = 3
                    ta = df_profile_max[tavar].values[0]

                    # WOCE quality control flags are used: 2 = good value, 3 = questionable value, 4 = bad value,
                    # 5 = value not reported, 6 = mean of replicate measurements, 9 = sample not drawn.
                    # use default TA value if TA is not good
                    if df_profile_max['TA_Flag'].values[0] not in (2, 6):
                        ta = 2200

                    par2_type = 1

                    kwargs = dict(salinity=df_profile_max.CTDSAL_PSS78.values[0],
                                  temperature=report_temp,
                                  temperature_out=df_profile_max.CTDTEMP_ITS90.values[0],
                                  pressure=0,
                                  pressure_out=df_profile_max.CTDPRES_dbar.values[0],
                                  opt_pH_scale=1,
                                  opt_k_carbonic=4,
                                  opt_k_bisulfate=1,
                                  opt_total_borate=1,
                                  opt_k_fluoride=2)

                    results = pyco2.sys(ph_20c, ta, par1_type, par2_type, **kwargs)
                    pH_bottom = np.round(results['pH_out'], 4)
                    temp_bottom = df_profile_max.CTDTEMP_ITS90.values[0]

                    # calculate aragonite
                    omega_bottom, pco2, revelle = cf.run_co2sys_ta_ph(df_profile_max[tavar].values[0],
                                                                      pH_bottom,
                                                                      df_profile_max.CTDSAL_PSS78.values[0],
                                                                      temp_bottom,
                                                                      df_profile_max.CTDPRES_dbar.values[0])

                    data['data_vars']['pH_bottom']['data'] = np.append(data['data_vars']['pH_bottom']['data'],
                                                                       pH_bottom)
                    data['data_vars']['omega_bottom']['data'] = np.append(data['data_vars']['omega_bottom']['data'],
                                                                          np.round(omega_bottom, 2))
                    data['data_vars']['temperature_bottom']['data'] = np.append(
                        data['data_vars']['temperature_bottom']['data'], temp_bottom)

    # save as netcdf
    outds = xr.Dataset.from_dict(data)

    # add created time to global attrs
    datetime_format = '%Y-%m-%dT%H:%M:%SZ'
    created = dt.datetime.now(dt.UTC).strftime(datetime_format)  # creation time Timestamp
    time_start = dt.datetime.fromtimestamp(np.nanmin(outds.time.values), dt.UTC).strftime('%Y-%m-%d')
    time_end = dt.datetime.fromtimestamp(np.nanmax(outds.time.values), dt.UTC).strftime('%Y-%m-%d')
    start_yr = dt.datetime.fromtimestamp(np.nanmin(outds.time.values), dt.UTC).strftime('%Y')
    end_yr = dt.datetime.fromtimestamp(np.nanmax(outds.time.values), dt.UTC).strftime('%Y')

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
        if k not in ['data_source', 'cruise', 'obs_type']:
            encoding[k] = {'zlib': True, 'complevel': 1, '_FillValue': -999}

    encoding['time'] = dict(zlib=False, _FillValue=False, dtype=np.double)

    save_file = os.path.join(savedir, f'vessel_based_bottom_OA_data_{start_yr}_{end_yr}.nc')
    outds.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims='time')

    # summarize the data sources
    ds = xr.open_dataset(save_file)
    ds['year'] = ds['time.year']
    df = ds.to_pandas()
    df.reset_index(inplace=True)
    summary = df.groupby(['data_source', 'cruise', 'obs_type', 'accession', 'year']).size()
    summary = summary.reset_index()
    summary.rename({0: 'sample_count'})
    summary.sort_values(['data_source', 'year'], inplace=True)
    save_file = os.path.join(savedir, f'vessel_based_bottom_OA_source_summary_{start_yr}_{end_yr}.csv')
    summary.to_csv(save_file, index=False)


if __name__ == '__main__':
    lons = [-78, -65, -65, -78]  # longitude boundaries for grabbing vessel-based data
    lats = [35, 35, 45, 45]  # latitude boundaries for grabbing vessel-based data
    codap = '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/CODAPv2021/CODAP_NA_v2021.nc'
    other = ['/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2019/Accession_0209156-discrete-profiles/33GG20190815-GU1902_data.csv',
              '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2019/Accession_0209045-discrete-profiles/33HH20190522-HB1902_data.csv',
              '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2019/Accession_0210238-discrete-profiles/33GG20191015-GU1905_data.csv',
              '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2021/Accession_0248269-discrete-profiles/33GG20210514-GU2102_data.csv',
              '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2021/Accession_0249432-discrete-profiles/334B20210805-PC2104_Data.csv',
              '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2021/Accession_0249517-discrete-profiles/334B20211015-PC2106_Data.csv',
              '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2022/Accession_0276023-discrete-profiles/33HH20220531_HB2204_Data.csv',
              '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2022/Accession_0283758-discrete-profiles/334B20221101_PC2205_Data-mod.csv',
              '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2023/Accession_0296717-discrete-profiles/33HH20230609_HB2302_Data.csv',
             '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/ECOA/ECOA-3/Accession_0283329-discrete/ECOA_3_CTD_MasterDataSheet_09_26_2023_Accession_0283329-mod.csv']
    save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/output_nc'
    main(lons, lats, codap, other, save_directory)
