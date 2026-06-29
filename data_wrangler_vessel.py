#!/usr/bin/env python

"""
Author: Lori Garzio on 10/18/2024
Last modified: 6/26/2026
Grab vessel-based surface- and bottom-water pH and omega data from CODAP-NA and additional ECOMON 
and ECOA datasets. Export as NetCDF.
CODAP-NA v2021 dataset documented here: https://essd.copernicus.org/articles/13/2777/2021/
Must run data_wrangler_codap_v2021.py first to create the CODAP-NA v2021 summary csv file.
Additional cruise datasets were downloaded from the NCEI OCADs data portal
(https://www.ncei.noaa.gov/products/ocean-carbon-acidification-data-system)
"""

import os
import yaml
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
    # read in global bathymetry file and subset to the region
    # for comparing profile bottom depth to actual bottom depth
    bathymetry = '/Users/lgarzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    extent = [-78, -65, 35, 45]
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))
    
    # get variable configuration file
    cfile = os.path.join(os.path.dirname(__file__), 'configs', 'vessel_vars.yml')
    with open(cfile, "r") as file:
        config_dict = yaml.safe_load(file)

    # get filepaths for other datasets not included in CODAP-NA v2021
    with open(extra_files, "r") as file:
        extra_files_list = yaml.safe_load(file)
    extra_files = extra_files_list.split(' ')

    with open(underway_files, "r") as file:
        underway_files_list = yaml.safe_load(file)
    underway_files = underway_files_list.split(' ')

    # initialize dictionary to append surface and bottom data from cruises (CODAPv2021, ECOMON and ECOA datasets)
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
            "comment": "Synthesis of surface and bottom pH and aragonite saturation state data from vessel-based measurements "
                       "that were spatially limited to the U.S. Northeast Shelf.",
            "data_sources": "CODAP-NA v2021 dataset (https://essd.copernicus.org/articles/13/2777/2021/) and subsequent "
                            "ECOMON and ECOA cruises not included in the CODAP-NA v2021 dataset. Additional datasets "
                            "were downloaded from NCEI's Ocean Carbon and Acidification Data Portal "
                            "(https://www.ncei.noaa.gov/access/ocean-carbon-acidification-data-system-portal/)"
        },
        "dims": "time",
        "data_vars": {}
    }
    
    # add variables from the config file to the data dictionary
    for varname, var_dict in config_dict.items():
        dtype = var_dict.pop('dtype')
        var_dict['data'] = np.array([], dtype=dtype)
        data['data_vars'][varname] = var_dict

#######################################################################################################################
    # get CODAP data
    df = pd.read_csv(codap_file)

    # for each profile on each cruise, find the shallowest and deepest samples
    cruises = np.unique(df.Cruise_ID)
    for cruise in cruises:
        dfc = df[df.Cruise_ID == cruise]
        profile_num = np.unique(dfc.Profile_number)
        for profile in profile_num:
            dfc_profile = dfc[dfc.Profile_number == profile]
            profile_coords = [dfc_profile.Longitude.values[0], dfc_profile.Latitude.values[0]]
            tm = pd.to_datetime(dfc_profile.time.values[0]).timestamp()
            
            # find surface sample
            mindepth, pH_surface, omega_surface, temp_surface = cf.return_surface_data(dfc_profile)

            # find deep sample
            maxdepth, pH_bottom, omega_bottom, temp_bottom, station_water_depth = cf.return_bottom_data(dfc_profile, bathy, profile_coords)


            # add data to dictionary
            data['coords']['time']['data'] = np.append(data['coords']['time']['data'], tm)

            data['data_vars']['data_source']['data'] = np.append(data['data_vars']['data_source']['data'],
                                                                    'CODAP_NA_v2021')
            data['data_vars']['cruise']['data'] = np.append(data['data_vars']['cruise']['data'], cruise)
            data['data_vars']['obs_type']['data'] = np.append(data['data_vars']['obs_type']['data'],
                                                                dfc_profile.Observation_type.values[0])
            data['data_vars']['accession']['data'] = np.append(data['data_vars']['accession']['data'],
                                                                int(dfc_profile.Accession.values[0]))
            
            data['data_vars']['depth_surface']['data'] = np.append(data['data_vars']['depth_surface']['data'], mindepth)
            data['data_vars']['depth_bottom']['data'] = np.append(data['data_vars']['depth_bottom']['data'], maxdepth)
            data['data_vars']['station_water_depth']['data'] = np.append(data['data_vars']['station_water_depth']['data'],
                                                                          station_water_depth)
            data['data_vars']['lat']['data'] = np.append(data['data_vars']['lat']['data'], profile_coords[1])
            data['data_vars']['lon']['data'] = np.append(data['data_vars']['lon']['data'], profile_coords[0])

            data['data_vars']['pH_surface']['data'] = np.append(data['data_vars']['pH_surface']['data'], pH_surface)
            data['data_vars']['omega_surface']['data'] = np.append(data['data_vars']['omega_surface']['data'], omega_surface)
            data['data_vars']['temperature_surface']['data'] = np.append(data['data_vars']['temperature_surface']['data'], temp_surface)
            data['data_vars']['pH_bottom']['data'] = np.append(data['data_vars']['pH_bottom']['data'], pH_bottom)
            data['data_vars']['omega_bottom']['data'] = np.append(data['data_vars']['omega_bottom']['data'], omega_bottom)
            data['data_vars']['temperature_bottom']['data'] = np.append(data['data_vars']['temperature_bottom']['data'], temp_bottom)

#######################################################################################################################
    # add extra datasets
    accession_mapping = {'HB1902': '0209045',
                         'GU1902': '0209156',
                         'GU1905': '0210238',
                         'GU2102': '0248269',
                         'PC2104': '0249432',
                         'PC2106': '0249517',
                         'PC2205': '0283758',
                         'HB2204': '0276023',
                         'ECOA3': '0283329',
                         'HB2302': '0296717',
                         'HB2303': '0302973',
                         'PC2305': '0303262',
                         'HB2401': '0309660',
                         'HB2403': '0309662',
                         'HB2406': '0309663',
                         'PC2406': '0309698'}

    # additional datasets that aren't included in CODAP v2021
    for ef in extra_files:
        print(ef) # TODO 33GG20210514-GU2102_data.csv should have flow-thru and niskin sampling ########
        df = pd.read_csv(ef)
        df.replace(-999, np.nan, inplace=True)

        # get cruise ID
        unique_cruise = np.unique(df.Cruise_ID)
        if len(unique_cruise) > 1:
            raise ValueError(f'More than one cruise ID found in {ef}')
        else:
            cruise = unique_cruise[0]
            cruise_accession = accession_mapping[cruise]

        try:
            df.dropna(subset=['pH_TS_20C'], inplace=True)
            phvar = 'pH_TS_20C'
            tavar = 'TA_umol/kg'
            df['report_temp'] = 20
        except KeyError:
            try:
                df.dropna(subset=['pH 20C'], inplace=True)
                phvar = 'pH 20C'
                tavar = 'TA (umol/kg)'
                df['report_temp'] = 20
            except KeyError:
                # ECOA-3 dataset
                # ECOA2 and ECOA3 - pH measurements were made at 25C (pers comm Shawn Shellito <Shawn.Shellito@unh.edu>)
                df.dropna(subset=['pH_T_measured_25C'], inplace=True)
                # get rid of questionable (3) and bad (4) pH flags
                df = df[df.pH_flag != 3]
                df = df[df.pH_flag != 4]

                phvar = 'pH_T_measured_25C'
                tavar = 'TA'
                df['report_temp'] = 25

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

        # calculate pH at in situ temperature, pressure, and salinity using PyCO2SYS
        pH_reported = df[phvar].values
        par1_type = 3

        # WOCE quality control flags are used: 2 = good value, 3 = questionable value, 4 = bad value,
        # 5 = value not reported, 6 = mean of replicate measurements, 9 = sample not drawn.

        # Use TA value if the TA flag is 2 (good) or 6 (mean of replicate measurements). 
        # Otherwise, use default TA value of 2200 if the TA flag is 3 (questionable), 4 (bad), 
        # 5 (not reported), or 9 (sample not drawn).
        ta = np.where(df['TA_Flag'].isin([2, 6]), df[tavar].values, 2200)
        par2_type = 1

        kwargs = dict(salinity=df.CTDSAL_PSS78.values,
                      temperature=df['report_temp'].values,
                      temperature_out=df.CTDTEMP_ITS90.values,
                      pressure=0,
                      pressure_out=df.CTDPRES_dbar.values,
                      opt_pH_scale=1,
                      opt_k_carbonic=4,
                      opt_k_bisulfate=1,
                      opt_total_borate=1,
                      opt_k_fluoride=2)

        results = pyco2.sys(pH_reported, ta, par1_type, par2_type, **kwargs)
        pH_insitu = np.round(results['pH_out'], 4)
        df['pH_insitu'] = pH_insitu

        # calculate aragonite
        omega, pco2, revelle = cf.run_co2sys_ta_ph(ta,
                                                   pH_insitu,
                                                   df.CTDSAL_PSS78.values,
                                                   df.CTDTEMP_ITS90.values,
                                                   df.CTDPRES_dbar.values)
        
        df['aragonite'] = np.round(omega, 2)

        # find the column name for sample depth
        sample_depth_vars = ['Depth_meters', 'Depth_sampling (M)', 'Depth']
        for sample_depth_var in sample_depth_vars:
            try:
                mindepth = np.nanmin(df[sample_depth_var])
                depthvar = sample_depth_var
            except KeyError:
                continue

        # find the column name for station water depth
        bottom_depth_vars = ['Depth_Bottom_meters', 'Depth_station (M)']
        bdv = None
        for bottom_depth_var in bottom_depth_vars:
            try:
                find_bottom_depth = np.nanmin(df[bottom_depth_var])
                bdv = bottom_depth_var
            except KeyError:
                continue

        # for each profile, find the shallowest sample
        profile_num = np.unique(df.profile_id)
        for profile in profile_num:
            df_profile = df[df.profile_id == profile]

            try:
                profile_coords = [df_profile.Longitude_Dec_Deg.values[0], df_profile.Latitude_Dec_Deg.values[0]]
            except AttributeError:
                profile_coords = [df_profile.Longitude.values[0], df_profile.Latitude.values[0]]
            tm = pd.to_datetime(df_profile.time.values[0]).timestamp()

            # find surface sample
            mindepth, pH_surface, omega_surface, temp_surface = cf.return_surface_data(df_profile,
                                                                                       depth_col=depthvar,
                                                                                       pH_col='pH_insitu',
                                                                                       omega_col='aragonite',
                                                                                       omega_est_col=None,
                                                                                       temp_col='CTDTEMP_ITS90')

            # find deep sample
            maxdepth, pH_bottom, omega_bottom, temp_bottom, station_water_depth = cf.return_bottom_data(df_profile, 
                                                                                                        bathy, 
                                                                                                        profile_coords,
                                                                                                        depth_col=depthvar,
                                                                                                        bottom_depth_col=bdv,
                                                                                                        pH_col='pH_insitu',
                                                                                                        omega_col='aragonite',
                                                                                                        omega_est_col=None,
                                                                                                        temp_col='CTDTEMP_ITS90')

            # add data to dictionary
            data['coords']['time']['data'] = np.append(data['coords']['time']['data'], tm)
            data['data_vars']['data_source']['data'] = np.append(data['data_vars']['data_source']['data'], 'NCEI-OCADS')
            data['data_vars']['cruise']['data'] = np.append(data['data_vars']['cruise']['data'], cruise)
            data['data_vars']['obs_type']['data'] = np.append(data['data_vars']['obs_type']['data'],
                                                                df_profile.Observation_Type.values[0])
            data['data_vars']['accession']['data'] = np.append(data['data_vars']['accession']['data'], cruise_accession)

            data['data_vars']['depth_surface']['data'] = np.append(data['data_vars']['depth_surface']['data'], mindepth)
            data['data_vars']['depth_bottom']['data'] = np.append(data['data_vars']['depth_bottom']['data'], maxdepth)
            data['data_vars']['station_water_depth']['data'] = np.append(data['data_vars']['station_water_depth']['data'],
                                                                          station_water_depth)
            data['data_vars']['lat']['data'] = np.append(data['data_vars']['lat']['data'], profile_coords[1])
            data['data_vars']['lon']['data'] = np.append(data['data_vars']['lon']['data'], profile_coords[0])

            data['data_vars']['pH_surface']['data'] = np.append(data['data_vars']['pH_surface']['data'], pH_surface)
            data['data_vars']['omega_surface']['data'] = np.append(data['data_vars']['omega_surface']['data'], omega_surface)
            data['data_vars']['temperature_surface']['data'] = np.append(data['data_vars']['temperature_surface']['data'], temp_surface)
            
            data['data_vars']['pH_bottom']['data'] = np.append(data['data_vars']['pH_bottom']['data'], pH_bottom)
            data['data_vars']['omega_bottom']['data'] = np.append(data['data_vars']['omega_bottom']['data'], omega_bottom)
            data['data_vars']['temperature_bottom']['data'] = np.append(data['data_vars']['temperature_bottom']['data'], temp_bottom)

########################################################################################################################
    # additional underway datasets
    accession_underway_mapping = {'ECOA-1': '0157389',
                                  'ECOA-2': '0215462',
                                  'ECOA3': '0295751'
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
        pH_reported = df[phvar].values
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

        results = pyco2.sys(pH_reported, ta, par1_type, par2_type, **kwargs)
        pH_insitu = np.round(results['pH_out'], 4)

        # calculate aragonite saturation state
        press = np.repeat(5, len(df.time))
        omega_surface, pco2, revelle = cf.run_co2sys_ta_ph(ta,
                                                           pH_insitu,
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

        data['data_vars']['depth_surface']['data'] = np.append(data['data_vars']['depth_surface']['data'], press)
        data['data_vars']['depth_bottom']['data'] = np.append(data['data_vars']['depth_bottom']['data'], 
                                                              np.repeat(np.nan, len(df.time)))
        data['data_vars']['station_water_depth']['data'] = np.append(data['data_vars']['station_water_depth']['data'],
                                                                        np.repeat(np.nan, len(df.time)))
        data['data_vars']['lat']['data'] = np.append(data['data_vars']['lat']['data'], np.array(df.Latitude))
        data['data_vars']['lon']['data'] = np.append(data['data_vars']['lon']['data'], np.array(df.Longitude))

        data['data_vars']['pH_surface']['data'] = np.append(data['data_vars']['pH_surface']['data'], pH_insitu)
        data['data_vars']['omega_surface']['data'] = np.append(data['data_vars']['omega_surface']['data'],
                                                               np.round(omega_surface, 2))
        data['data_vars']['temperature_surface']['data'] = np.append(data['data_vars']['temperature_surface']['data'], 
                                                                     df[tempvar].values)
        
        data['data_vars']['pH_bottom']['data'] = np.append(data['data_vars']['pH_bottom']['data'],
                                                           np.repeat(np.nan, len(df.time)))
        data['data_vars']['omega_bottom']['data'] = np.append(data['data_vars']['omega_bottom']['data'],
                                                              np.repeat(np.nan, len(df.time)))
        data['data_vars']['temperature_bottom']['data'] = np.append(data['data_vars']['temperature_bottom']['data'],
                                                                    np.repeat(np.nan, len(df.time)))

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
        if k not in ['data_source', 'cruise', 'obs_type', 'accession']:
            encoding[k] = {'zlib': True, 'complevel': 1, '_FillValue': -999}

    encoding['time'] = dict(zlib=False, _FillValue=False, dtype=np.double)

    save_file = os.path.join(savedir, f'TEST_vessel_based_OA_data_{start_yr}_{end_yr}.nc')
    outds.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims='time')

    # summarize the data sources
    ds = xr.open_dataset(save_file)
    ds['year'] = ds['time.year']
    df = ds.to_pandas()
    df.reset_index(inplace=True)
    summary = df.groupby(['data_source', 'cruise', 'obs_type', 'accession', 'year'], as_index=False).size()
    summary.sort_values(['data_source', 'year'], inplace=True)
    save_file = os.path.join(savedir, f'TEST_vessel_based_OA_source_summary_{start_yr}_{end_yr}.csv')
    summary.to_csv(save_file, index=False)


if __name__ == '__main__':
    lons = [-78, -65, -65, -78]  # longitude boundaries for grabbing vessel-based data
    lats = [35, 35, 45, 45]  # latitude boundaries for grabbing vessel-based data
    codap = '/Users/lgarzio/Documents/rucool/Saba/OA_cruise_data/CODAPv2021/CODAP-NA_v2021-filtered-for-analysis.csv'
    other = '/Users/lgarzio/Documents/rucool/Saba/OA_cruise_data/soe_other_filepaths.yml'
    other_underway = '/Users/lgarzio/Documents/rucool/Saba/OA_cruise_data/soe_other_underway_filepaths.yml'
    save_directory = '/Users/lgarzio/Documents/rucool/Saba/NOAA_SOE/data/output_nc'
    main(lons, lats, codap, other, other_underway, save_directory)
