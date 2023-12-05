#!/usr/bin/env python

"""
Author: Lori Garzio on 10/18/2021
Last modified: 11/17/2023
Grab summer bottom-water pH and omega data from glider datasets and export as NetCDF.
CODAP-NA dataset documented here: https://essd.copernicus.org/articles/13/2777/2021/
"""

import os
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
import statistics
from collections import OrderedDict
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console


def main(glider_files, savedir):
    bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    extent = [-78, -65, 35, 45]
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))

    # initialize dictionary to append bottom pH and aragonite data from glider deployments
    data = {
        "coords": {
            "time": {"dims": "time", "data": np.array([], dtype='datetime64[ns]')}
        },
        "attrs": {
            "comment": "Synthesis of bottom pH and aragonite saturation state data from glider-based measurements "
                       "collected during summer months (June-August) that were spatially limited to the U.S. Northeast "
                       "Shelf. "
        },
        "dims": "time",
        "data_vars": {
            "deployment": {
                "dims": "time",
                "data": np.array([], dtype='<U32'),
                "attrs": {
                    "units": "1",
                    "comment": "Source of data"
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
                    "comment": "Depth at which sample was collected"
                }
            },
            "bottom_depth": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "m",
                    "long_name": "Bottom Depth",
                    "description": "Depth of water column at the sampling location",
                    "comment": "Bottom depth was calculated using the sample coordinates and GEBCO’s gridded global "
                               "bathymetry (https://download.gebco.net/)"
                }
            },
            "pH_bottom": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "1",
                    "long_name": "pH",
                    "description": "Bottom pH",
                    "comment": "Bottom pH is the the median of the values within the deepest 1m of a glider profile, "
                               "provided the sampling depth was no shallower than the bottom 20% of total "
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
                    "comment": "Calculated using PyCO2SYS (Humphreys et al. 2020) with inputs of "
                               "pressure, temperature, salinity, total alkalinity, and pH. "
                               "Bottom aragonite saturation state is the the median of the values within the deepest "
                               "1m of a glider profile, provided the sampling depth was no shallower than the bottom "
                               "20% of total water column depth."
                }
            },
        },
    }

    # get glider data
    # for each glider, grab the data in the bottom 1m of each profile except when the profile is deeper than a
    # defined threshold (when the glider is off the shelf and no longer sampling the bottom)

    max_depth_threshold = dict(sbu01=170,
                               ru30=170,
                               um_242=300,
                               ru39=170)

    for gf in glider_files:
        ds = xr.open_dataset(gf)
        deployment = ds.title.split(' ')[0].split('-delayed')[0]
        glider = (deployment.split('-')[0]).lower()

        # grab the correct variables for each dataset
        if deployment in ['SBU01-20220805T1855', 'SBU01-20230706T1707']:  # datasets from the DAC
            profileid = np.unique(ds.time.values)
            profile_identifier = 'time'
            phvar = 'pHtot'
            arag = ds.aragonite_saturation_state.values
        else:
            try:
                profileid = np.unique(ds.profile_time.values)
                profile_identifier = 'profile_time'
                try:
                    phvar = 'ph_total_shifted'
                    ph = ds[phvar].values
                except KeyError:
                    phvar = 'pH_shifted'
            except AttributeError:
                profileid = np.unique(ds.profile_id.values)
                profile_identifier = 'profile_id'
                phvar = 'pH'
            try:
                arag = ds.saturation_aragonite.values
            except AttributeError:
                arag = ds.aragonite_saturation_state.values
        for pid in profileid:
            pidx = np.where(ds[profile_identifier].values == pid)[0]
            max_profile_depth = np.nanmax(ds.depth.values[pidx])
            if np.logical_and(max_profile_depth > 10, max_profile_depth < max_depth_threshold[glider]):
                idx = np.where(ds.depth.values[pidx] > max_profile_depth - 1)[0]

                # compare the glider depth to the global bathymetry file
                gl_lat = statistics.median(ds.latitude.values[pidx][idx])
                gl_lon = statistics.median(ds.longitude.values[pidx][idx])
                profile_coords = [gl_lon, gl_lat]
                lat_idx = abs(bathy.lat.values - profile_coords[1]).argmin()
                lon_idx = abs(bathy.lon.values - profile_coords[0]).argmin()
                station_water_depth = -bathy.elevation[lat_idx, lon_idx].values

                # if the glider sample is within +/- 20% of the water column, keep the value
                depth_threshold = [station_water_depth * .8, station_water_depth * 1.2]
                if np.logical_and(max_profile_depth > depth_threshold[0], max_profile_depth < depth_threshold[1]):
                    if np.sum(~np.isnan(ds[phvar].values[pidx][idx])) > 0:
                        tm = ds.time.values[pidx][idx][0]
                        pH_bottom = np.round(statistics.median(ds[phvar].values[pidx][idx]), 4)
                        omega_bottom = np.round(statistics.median(arag[pidx][idx]), 2)
                        data['coords']['time']['data'] = np.append(data['coords']['time']['data'], tm)
                        data['data_vars']['deployment']['data'] = np.append(data['data_vars']['deployment']['data'],
                                                                            deployment)
                        data['data_vars']['depth']['data'] = np.append(data['data_vars']['depth']['data'],
                                                                       statistics.median(ds.depth.values[pidx][idx]))
                        data['data_vars']['bottom_depth']['data'] = np.append(data['data_vars']['bottom_depth']['data'],
                                                                              station_water_depth)
                        data['data_vars']['lat']['data'] = np.append(data['data_vars']['lat']['data'], gl_lat)
                        data['data_vars']['lon']['data'] = np.append(data['data_vars']['lon']['data'], gl_lon)
                        data['data_vars']['pH_bottom']['data'] = np.append(data['data_vars']['pH_bottom']['data'],
                                                                           pH_bottom)
                        data['data_vars']['omega_bottom']['data'] = np.append(data['data_vars']['omega_bottom']['data'],
                                                                              omega_bottom)

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

    save_file = os.path.join(savedir, 'glider_based_summer_bottom_OA_data.nc')
    outds.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims='time')


if __name__ == '__main__':
    gliders = [
        '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/nc_files/sbu01-20210720T1628-profile-sci-delayed-qc_shifted_co2sys_final.nc',
        '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/nc_files/ru30-20190717T1812-profile-sci-delayed-dac.nc',
        '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/nc_files/ru30-20210716T1804-profile-sci-delayed_qc_final.nc',
        '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/nc_files/um_242-20210630T1916-profile-sci-delayed_shifted_final.nc',
        '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/nc_files/SBU01-20220805T1855-delayed-co2sys.nc',
        '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/nc_files/SBU01-20230706T1707-delayed-co2sys.nc',
        '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/nc_files/ru39-20230817T1520-delayed.nc'
    ]
    save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/output_nc'
    main(gliders, save_directory)
