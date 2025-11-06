#!/usr/bin/env python

"""
Author: Lori Garzio on 10/18/2024
Last modified: 11/6/2025
Grab vessel-based surface-water pH and omega data from CODAP-NA and additional ECOMON and ECOA datasets. Export as NetCDF.
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


def main(lon_bounds, lat_bounds, codap_file, extra_files, underway_files, savedir):
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
            "comment": "Synthesis of surface pH and aragonite saturation state data from vessel-based measurements "
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
            "pH_surface": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "1",
                    "long_name": "Surface pH",
                    "comment": "For the CODAP-NA v2021 dataset, pH_TS_insitu_measured was used when available. "
                               "Otherwise, pH_TS_insitu_calculated was used. For the datasets downloaded from "
                               "NCEI, pH at in situ temperature, pressure, and salinity were calculated from "
                               "pH_TS_20C using PyCO2SYS.",
                    "description": "Surface pH is the shallowest measurement of a vertical CTD/Rosette cast "
                                   "if the sample depth was < 10m. For stations with multiple samples at the "
                                   "shallowest depth, the median value was calculated.",
                }
            },
            "omega_surface": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "1",
                    "long_name": "Surface Aragonite Saturation State",
                    "description": "Surface aragonite saturation state is the shallowest measurement of a vertical "
                                   "CTD/Rosette cast if the sample depth was < 10m. For stations with multiple "
                                   "samples at the shallowest depth, the median value was calculated.",
                    "comment": "Aragonite saturation state from the original data source, when available. "
                               "When not available, calculated using PyCO2SYS (Humphreys et al. 2020) with inputs of "
                               "pressure, temperature, salinity, total alkalinity, and pH. "

                }
            },
            "temperature_surface": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "C",
                    "long_name": "Surface Temperature",
                    "description": "Surface seawater temperature is the shallowest measurement of a vertical "
                                   "CTD/Rosette cast if the sample depth was < 10m. For stations with multiple "
                                   "samples at the shallowest depth, the median value was calculated."
                }
            },
        },
    }

#######################################################################################################################
    # get CODAP data
    try:
        # see if the summary csv file already exists
        df = pd.read_csv(os.path.join(os.path.dirname(codap_file), 'CODAP-NA_v2021-filtered-for-surface-analysis.csv'))
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
        df.to_csv(os.path.join(os.path.dirname(codap_file), 'CODAP-NA_v2021-filtered-for-surface-analysis.csv'))

    # for each profile on each cruise, find the shallowest sample
    cruises = np.unique(df.Cruise_ID)
    for cruise in cruises:
        dfc = df[df.Cruise_ID == cruise]
        profile_num = np.unique(dfc.Profile_number)
        for profile in profile_num:
            dfc_profile = dfc[dfc.Profile_number == profile]
            mindepth = np.nanmin(dfc_profile.Depth)
            if bool(mindepth < 10):
                dfc_profile_min = dfc_profile[dfc_profile.Depth == mindepth]
                profile_coords = [dfc_profile_min.Longitude.values[0], dfc_profile_min.Latitude.values[0]]

                # drop lines where measured omega isn't available
                if len(dfc_profile_min) > 1:
                    print('check')
                    dfc_profile_min = dfc_profile_min[dfc_profile_min.Aragonite > 0]

                    # if you removed all rows of data, go back to the original
                    if len(dfc_profile_min) < 1:
                        dfc_profile_min = dfc_profile[dfc_profile.Depth == mindepth]

                tm = pd.to_datetime(dfc_profile_min.time.values[0]).timestamp()
                pH_surface = np.nanmedian(np.array(dfc_profile_min.pH))

                omega_surface = np.nanmedian(np.array(dfc_profile_min.Aragonite))

                # if measured omega isn't available (-999) use estimated aragonite
                if bool(omega_surface < 0):
                    omega_surface = np.nanmedian(np.array(dfc_profile_min.aragonite_estimated))

                temp_surface = np.nanmedian(np.array(dfc_profile_min.CTDTEMP_ITS90))
                if bool(temp_surface < 0):
                    temp_surface = np.nan

                # add data to dictionary
                data['coords']['time']['data'] = np.append(data['coords']['time']['data'], tm)

                data['data_vars']['data_source']['data'] = np.append(data['data_vars']['data_source']['data'],
                                                                     'CODAP_NA_v2021')
                data['data_vars']['cruise']['data'] = np.append(data['data_vars']['cruise']['data'], cruise)
                data['data_vars']['obs_type']['data'] = np.append(data['data_vars']['obs_type']['data'],
                                                                  dfc_profile_min.Observation_type.values[0])
                data['data_vars']['accession']['data'] = np.append(data['data_vars']['accession']['data'],
                                                                   int(dfc_profile_min.Accession.values[0]))
                data['data_vars']['depth']['data'] = np.append(data['data_vars']['depth']['data'], mindepth)
                data['data_vars']['lat']['data'] = np.append(data['data_vars']['lat']['data'], profile_coords[1])
                data['data_vars']['lon']['data'] = np.append(data['data_vars']['lon']['data'], profile_coords[0])
                data['data_vars']['pH_surface']['data'] = np.append(data['data_vars']['pH_surface']['data'], pH_surface)
                data['data_vars']['omega_surface']['data'] = np.append(data['data_vars']['omega_surface']['data'],
                                                                       omega_surface)
                data['data_vars']['temperature_surface']['data'] = np.append(
                    data['data_vars']['temperature_surface']['data'], temp_surface)

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
                         'HB2302': 296717,
                         'HB2303': 302973,
                         'PC2305': 303262}

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

        # for each profile, find the shallowest sample
        profile_num = np.unique(df.profile_id)
        for profile in profile_num:
            df_profile = df[df.profile_id == profile]

            # find the column name for sample depth
            sample_depth_vars = ['Depth_meters', 'Depth_sampling (M)', 'Depth']
            for sample_depth_var in sample_depth_vars:
                try:
                    mindepth = np.nanmin(df_profile[sample_depth_var])
                    depthvar = sample_depth_var
                except KeyError:
                    continue

            # depth must be <10m
            if bool(mindepth < 10):

                df_profile_min = df_profile[df_profile[depthvar] == mindepth]

                try:
                    profile_coords = [df_profile_min.Longitude_Dec_Deg.values[0],
                                      df_profile_min.Latitude_Dec_Deg.values[0]]
                except AttributeError:
                    profile_coords = [df_profile_min.Longitude.values[0],
                                      df_profile_min.Latitude.values[0]]

                tm = pd.to_datetime(df_profile_min['time'].values[0]).timestamp()
                cruise = df_profile_min.Cruise_ID.values[0]

                data['coords']['time']['data'] = np.append(data['coords']['time']['data'], tm)
                data['data_vars']['data_source']['data'] = np.append(data['data_vars']['data_source']['data'],
                                                                     'NCEI-OCADS')
                data['data_vars']['cruise']['data'] = np.append(data['data_vars']['cruise']['data'], cruise)
                data['data_vars']['obs_type']['data'] = np.append(data['data_vars']['obs_type']['data'],
                                                                  df_profile_min.Observation_Type.values[0])
                cruise_accession = accession_mapping[cruise]
                data['data_vars']['accession']['data'] = np.append(data['data_vars']['accession']['data'], cruise_accession)

                data['data_vars']['depth']['data'] = np.append(data['data_vars']['depth']['data'], mindepth)
                data['data_vars']['lat']['data'] = np.append(data['data_vars']['lat']['data'], profile_coords[1])
                data['data_vars']['lon']['data'] = np.append(data['data_vars']['lon']['data'], profile_coords[0])

                # calculate pH at in situ temperature, pressure, and salinity using PyCO2SYS
                ph_20c = df_profile_min[phvar].values[0]
                par1_type = 3
                ta = df_profile_min[tavar].values[0]

                # WOCE quality control flags are used: 2 = good value, 3 = questionable value, 4 = bad value,
                # 5 = value not reported, 6 = mean of replicate measurements, 9 = sample not drawn.
                # use default TA value if TA is not good
                if df_profile_min['TA_Flag'].values[0] not in (2, 6):
                    ta = 2200

                par2_type = 1

                kwargs = dict(salinity=df_profile_min.CTDSAL_PSS78.values[0],
                              temperature=report_temp,
                              temperature_out=df_profile_min.CTDTEMP_ITS90.values[0],
                              pressure=0,
                              pressure_out=df_profile_min.CTDPRES_dbar.values[0],
                              opt_pH_scale=1,
                              opt_k_carbonic=4,
                              opt_k_bisulfate=1,
                              opt_total_borate=1,
                              opt_k_fluoride=2)

                results = pyco2.sys(ph_20c, ta, par1_type, par2_type, **kwargs)
                pH_surface = np.round(results['pH_out'], 4)
                temp_surface = df_profile_min.CTDTEMP_ITS90.values[0]

                # calculate aragonite
                omega_surface, pco2, revelle = cf.run_co2sys_ta_ph(df_profile_min[tavar].values[0],
                                                                  pH_surface,
                                                                  df_profile_min.CTDSAL_PSS78.values[0],
                                                                  temp_surface,
                                                                  df_profile_min.CTDPRES_dbar.values[0])

                data['data_vars']['pH_surface']['data'] = np.append(data['data_vars']['pH_surface']['data'],
                                                                   pH_surface)
                data['data_vars']['omega_surface']['data'] = np.append(data['data_vars']['omega_surface']['data'],
                                                                      np.round(omega_surface, 2))
                data['data_vars']['temperature_surface']['data'] = np.append(
                    data['data_vars']['temperature_surface']['data'], temp_surface)

########################################################################################################################
    # additional underway datasets
    accession_underway_mapping = {'ECOA-1': 157389,
                                  'ECOA-2': 215462,
                                  'ECOA3': 295751
                                  }

    for uf in underway_files:
        print(uf)
        df = pd.read_csv(uf)
        df.replace(-999, np.nan, inplace=True)
        try:
            df.dropna(subset=['PH_TOT'], inplace=True)
            phvar = 'PH_TOT'
            phflagvar = 'PH_TOT_FLAG_W'
            tavar = 'ALKALI_umol_kg'
            taflagvar = 'ALKALI_FLAG_W'
            salvar = 'TSG_Sal'
            tempvar = 'TSG_SST_C'
        except KeyError:
            # ECOA-3 dataset
            df.dropna(subset=['pH_T_measured'], inplace=True)
            phvar = 'pH_T_measured'
            phflagvar = 'pH_flag'
            tavar = 'TA_umol_kg'
            taflagvar = 'TA_flag'
            salvar = 'SSS'
            tempvar = 'TSG_temp_ITS90'

        # WOCE quality control flags are used: 2 = good value, 3 = questionable value, 4 = bad value,
        # 5 = value not reported, 6 = mean of replicate measurements, 9 = sample not drawn.
        # use default TA value if TA is not good

        # remove questionable (3) and bad (4) pH
        df = df[df[phflagvar] != 3]
        df = df[df[phflagvar] != 4]

        # remove questionable (3) and bad (4) TA
        df = df[df[taflagvar] != 3]
        df = df[df[taflagvar] != 4]

        # format date
        try:
            df['time'] = pd.to_datetime(df[['Year', 'Month', 'Day']])
        except KeyError:
            df['time'] = df['Date_UTC'].map(lambda t: pd.to_datetime(str(t)))

        # make sure the data are within the defined extent
        df['in_region'] = ''
        for i, row in df.iterrows():
            if Polygon(list(zip(lon_bounds, lat_bounds))).contains(Point(row.Longitude, row.Latitude)):
                df.loc[i, 'in_region'] = 'yes'
            else:
                df.loc[i, 'in_region'] = 'no'

        # drop data if it's not in the region specified
        df = df[df.in_region == 'yes']

        # calculate pH at in situ temperature, pressure, and salinity using PyCO2SYS
        # assume pH is reported at 20C unless otherwise specified
        ph_20c = df[phvar].values
        par1_type = 3
        ta = df[tavar].values
        par2_type = 1

        # ECOA2 and ECOA3 - pH measurements were made at 25C (pers comm Shawn Shellito <Shawn.Shellito@unh.edu>)
        reported_temp = df.PH_TOT_TEMP_C.values

        kwargs = dict(salinity=df[salvar].values,
                      temperature=reported_temp,  # reported measurement temperature
                      temperature_out=df[tempvar].values,
                      pressure=0,
                      pressure_out=5,
                      opt_pH_scale=1,
                      opt_k_carbonic=4,
                      opt_k_bisulfate=1,
                      opt_total_borate=1,
                      opt_k_fluoride=2)

        results = pyco2.sys(ph_20c, ta, par1_type, par2_type, **kwargs)
        pH_surface = np.round(results['pH_out'], 4)

        # calculate aragonite saturation state
        press = np.repeat(5, len(df.time))
        omega_surface, pco2, revelle = cf.run_co2sys_ta_ph(ta,
                                                           pH_surface,
                                                           df[salvar].values,
                                                           df[tempvar].values,
                                                           press)

        # convert time to timestamp
        tm = df['time'].map(lambda t: pd.to_datetime(t).timestamp())

        # get cruise
        try:
            cruise = np.array(df.Cruise_ID)
        except AttributeError:
            cruise = np.array(df.CRUISE_ID)

        data['coords']['time']['data'] = np.append(data['coords']['time']['data'], np.array(tm))
        data['data_vars']['data_source']['data'] = np.append(data['data_vars']['data_source']['data'],
                                                             np.repeat('NCEI-OCADS', len(df.time)))
        data['data_vars']['cruise']['data'] = np.append(data['data_vars']['cruise']['data'], cruise)
        data['data_vars']['obs_type']['data'] = np.append(data['data_vars']['obs_type']['data'],
                                                          np.repeat('Flow-through', len(df.time)))
        cruise_accession = accession_underway_mapping[np.unique(cruise)[0]]
        data['data_vars']['accession']['data'] = np.append(data['data_vars']['accession']['data'],
                                                           np.repeat(cruise_accession, len(df.time)))

        data['data_vars']['depth']['data'] = np.append(data['data_vars']['depth']['data'], press)
        data['data_vars']['lat']['data'] = np.append(data['data_vars']['lat']['data'], np.array(df.Latitude))
        data['data_vars']['lon']['data'] = np.append(data['data_vars']['lon']['data'], np.array(df.Longitude))

        data['data_vars']['pH_surface']['data'] = np.append(data['data_vars']['pH_surface']['data'],
                                                            pH_surface)
        data['data_vars']['omega_surface']['data'] = np.append(data['data_vars']['omega_surface']['data'],
                                                               np.round(omega_surface, 2))
        data['data_vars']['temperature_surface']['data'] = np.append(
            data['data_vars']['temperature_surface']['data'], df[tempvar].values)

    # save as netcdf
    outds = xr.Dataset.from_dict(data)

    # add created time to global attrs
    datetime_format = '%Y-%m-%dT%H:%M:%SZ'
    created = dt.datetime.now(dt.UTC).strftime(datetime_format) # creation time Timestamp
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

    save_file = os.path.join(savedir, f'vessel_based_surface_OA_data_{start_yr}_{end_yr}.nc')
    outds.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims='time')

    # summarize the data sources
    ds = xr.open_dataset(save_file)
    ds['year'] = ds['time.year']
    df = ds.to_pandas()
    df.reset_index(inplace=True)
    summary = df.groupby(['data_source', 'cruise', 'obs_type', 'accession', 'year'], as_index=False).size()
    summary.sort_values(['data_source', 'year'], inplace=True)
    save_file = os.path.join(savedir, f'vessel_based_surface_OA_source_summary_{start_yr}_{end_yr}.csv')
    summary.to_csv(save_file, index=False)


if __name__ == '__main__':
    lons = [-78, -65, -65, -78]  # longitude boundaries for grabbing vessel-based data
    lats = [35, 35, 45, 45]  # latitude boundaries for grabbing vessel-based data
    codap = '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/CODAPv2021/CODAP_NA_v2021.nc'
    other = [
        '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2019/Accession_0209156-discrete-profiles/33GG20190815-GU1902_data.csv',
        '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2019/Accession_0209045-discrete-profiles/33HH20190522-HB1902_data.csv',
        '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2019/Accession_0210238-discrete-profiles/33GG20191015-GU1905_data.csv',
        '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2021/Accession_0248269-discrete-profiles/33GG20210514-GU2102_data.csv',
        '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2021/Accession_0249432-discrete-profiles/334B20210805-PC2104_Data.csv',
        '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2021/Accession_0249517-discrete-profiles/334B20211015-PC2106_Data.csv',
        '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2022/Accession_0276023-discrete-profiles/33HH20220531_HB2204_Data.csv',
        '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2022/Accession_0283758-discrete-profiles/334B20221101_PC2205_Data-mod.csv',
        '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2023/Accession_0296717-discrete-profiles/33HH20230609_HB2302_Data.csv',
        '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/ECOA/ECOA-3/Accession_0283329-discrete/ECOA_3_CTD_MasterDataSheet_09_26_2023_Accession_0283329-mod.csv',
        '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2023/Accession_0302973-discrete-profiles/33HH20230808_HB2303_Data.csv',
        '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/EcoMon/2023/Accession_0303262-discrete-profiles/334B20231027_PC2305_Data.csv'
        ]
    other_underway = [
        '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/ECOA/ECOA-1/Accession_0157389-underway/Discrete_Underway_Data_12082016_Accession_0157389-mod.csv',
        '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/ECOA/ECOA-2/Accession_0196423-surface-discrete/ECOA2_Discrete_Underway_Data-mod.csv',
        '/Users/garzio/Documents/rucool/Saba/OA_cruise_data/ECOA/ECOA-3/Accession_0295751-underway/ECOA_3_Underway_MasterDataSheet_NCEI_20240731._Accession_0295751-mod.csv'
        ]
    save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/output_nc'
    main(lons, lats, codap, other, other_underway, save_directory)
