#!/usr/bin/env python

"""
Author: Lori Garzio on 10/23/2024
Last modified: 10/23/2024
Grab bottom- and surface-water pH and omega data from glider datasets and export as NetCDF.
Datasets are available on the IOOS Glider DAC ERDDAP server https://gliders.ioos.us/erddap/index.html
"""

import os
import glob
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
import statistics
from collections import OrderedDict
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console
np.set_printoptions(suppress=True)


def interpolate_var(da):
    df = da.to_dataframe()
    var_interp = df[da.name].interpolate(method='linear', limit_direction='both', limit=2).values
    return var_interp


def main(rufiles, rufiles2, sbufiles, umfiles, savedir):
    bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    extent = [-78, -65, 35, 45]
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))

    # initialize dictionary to append bottom pH and aragonite data from glider deployments
    bottom_data = {
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
            "comment": "Synthesis of bottom pH and aragonite saturation state data from glider-based measurements "
                       "that were spatially limited to the U.S. Northeast Shelf."
        },
        "dims": "time",
        "data_vars": {
            "deployment": {
                "dims": "time",
                "data": np.array([], dtype='object'),
                "attrs": {
                    "units": "1",
                    "comment": "Deployment name"
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
                    "long_name": "Depth"
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
                    "long_name": "Bottom pH",
                    "description": "Bottom pH is the the median of the values within the deepest 1m of a glider "
                                   "profile, provided the sampling depth was no shallower than the bottom 20% of total "
                                   "water column depth."
                }
            },
            "omega_bottom": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "1",
                    "long_name": "Bottom Aragonite Saturation State",
                    "description": "Bottom aragonite saturation state is the the median of the values within the deepest "
                               "1m of a glider profile, provided the sampling depth was no shallower than the bottom "
                               "20% of total water column depth.",
                    "comment": "Calculated using PyCO2SYS (Humphreys et al. 2020) with inputs of "
                               "pressure, temperature, salinity, total alkalinity, and pH."
                }
            },
            "temperature_bottom": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "C",
                    "long_name": "Bottom Temperature",
                    "description": "Bottom temperature is the the median of the values within the deepest "
                                   "1m of a glider profile, provided the sampling depth was no shallower than the bottom "
                                   "20% of total water column depth."
                }
            },
        },
    }

    # initialize dictionary to append surface pH and aragonite data from glider deployments
    surface_data = {
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
            "comment": "Synthesis of surface pH and aragonite saturation state data from glider-based measurements "
                       "that were spatially limited to the U.S. Northeast Shelf."
        },
        "dims": "time",
        "data_vars": {
            "deployment": {
                "dims": "time",
                "data": np.array([], dtype='object'),
                "attrs": {
                    "units": "1",
                    "comment": "Deployment name"
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
                    "long_name": "Depth"
                }
            },
            "pH_surface": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "1",
                    "long_name": "Surface pH",
                    "description": "Surface pH is the the median of the values recorded at the top of a glider profile "
                                   "(between 2-4m depth), provided the glider profile reached at least 10m depth."
                }
            },
            "omega_surface": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "1",
                    "long_name": "Surface Aragonite Saturation State",
                    "description": "Surface aragonite saturation state is the the median of the values recorded at the "
                                   "top of a glider profile (between 2-4m depth), provided the glider profile reached at least 10m depth.",
                    "comment": "Calculated using PyCO2SYS (Humphreys et al. 2020) with inputs of "
                               "pressure, temperature, salinity, total alkalinity, and pH."
                }
            },
            "temperature_surface": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "C",
                    "long_name": "Surface Temperature",
                    "description": "Surface temperature is the the median of the values recorded at the "
                                   "top of a glider profile (between 2-4m depth), provided the glider profile reached at least 10m depth."
                }
            },
        },
    }

    # get glider data
    glider_files = []
    for f in glob.glob(os.path.join(rufiles, '*.nc')):
        glider_files.append(f)
    for f in glob.glob(os.path.join(rufiles2, '*.nc')):
        glider_files.append(f)
    for f in glob.glob(os.path.join(sbufiles, '*.nc')):
        glider_files.append(f)
    for f in glob.glob(os.path.join(umfiles, '*.nc')):
        glider_files.append(f)

    for gf in glider_files:
        print(gf)
        ds = xr.open_dataset(gf)
        deployment = ds.title.split(' ')[0].split('-delayed')[0]
        glider = (deployment.split('-')[0]).lower()

        if deployment == 'um_242-20210630T1904':
            ds = ds.swap_dims({'row': 'time'})
            ds = ds.sel(time=slice('2021-07-01', '2021-07-22'))

        # grab the correct variables for each dataset
        try:
            profile_identifier = 'profile_id'
            profileid = np.unique(ds[profile_identifier].values)
        except KeyError:
            profile_identifier = 'profile_time'
            profileid = np.unique(ds[profile_identifier].values)

        try:
            phvar = 'pH'
            ph = ds[phvar]
        except KeyError:
            try:
                phvar = 'pH_corrected'
                ph = ds[phvar]
            except KeyError:
                phvar = 'pHtot'
                ph = ds[phvar]

        try:
            latvar = 'latitude'
            lonvar = 'longitude'
            lat = ds[latvar]
            lon = ds[lonvar]
        except KeyError:
            latvar = 'profile_lat'
            lonvar = 'profile_lon'
            lat = ds[latvar]
            lon = ds[lonvar]

        # apply pH QC variables, if available
        qcvars = [f'{phvar}_qartod_gross_range_test', f'{phvar}_qartod_spike_test']
        for qv in qcvars:
            try:
                qc_idx = np.where(np.logical_or(ds[qv].values == 3, ds[qv].values == 4))[0]
                if len(qc_idx) > 0:
                    ds[phvar][qc_idx] = np.nan
            except KeyError:
                continue

        # calculate TA and omega for sbu gliders
        if glider in ['sbu01', 'sbu02']:
            sal_interp = interpolate_var(ds.salinity)
            temp_interp = interpolate_var(ds.temperature)
            press_interp = interpolate_var(ds.pressure)
            season = str(ds['time.season'].values[0])
            ta = cf.calc_ta_nyb(season, sal_interp)

            # when pH is nan, TA is nan
            idx = np.isnan(ds[phvar].values)
            ta[idx] = np.nan

            # run CO2SYS to calculate omega
            arag, pco2, revelle = cf.run_co2sys_ta_ph(ta,
                                                      ds[phvar].values,
                                                      sal_interp,
                                                      temp_interp,
                                                      press_interp)

        else:
            arag = ds.aragonite_saturation_state.values

        # apply QC to umaine gliders for some issues I noticed
        if glider in ['um_242']:
            # TA
            idx = np.where(ds.total_alkalinity.values < 1900)[0]
            ds[phvar].values[idx] = np.nan
            arag[idx] = np.nan

            # omega
            idx = np.where(arag > 5)[0]
            ds[phvar].values[idx] = np.nan
            arag[idx] = np.nan

            # pH
            idx = np.where(ds[phvar] > 8.6)[0]
            ds[phvar].values[idx] = np.nan
            arag[idx] = np.nan

        try:
            depth_interp = ds.depth_interpolated.values
        except AttributeError:
            depth_interp = interpolate_var(ds.depth)

        # for each glider, grab the data in the bottom 1m of each profile when the glider is within +/- 20% of the water
        # column depth (determined by comparison to global bathymetry data)
        # also grab the surface data (2-4m depth)
        for pid in profileid:
            pidx = np.where(ds[profile_identifier].values == pid)[0]
            max_profile_depth = np.nanmax(ds.depth.values[pidx])
            if max_profile_depth > 10:
                # grab surface data
                idx = np.where(np.logical_and(depth_interp[pidx] >= 2, depth_interp[pidx] <= 4))[0]
                if np.sum(~np.isnan(ds[phvar].values[pidx][idx])) > 0:
                    tm = pd.to_datetime(ds.time.values[pidx][idx][0]).timestamp()
                    pH_surface = np.round(np.nanmedian(ds[phvar].values[pidx][idx]), 4)
                    omega_surface = np.round(np.nanmedian(arag[pidx][idx]), 2)
                    temperature_surface = np.round(np.nanmedian(ds['temperature'].values[pidx][idx]), 2)
                    gl_lat = statistics.median(ds[latvar].values[pidx][idx])
                    gl_lon = statistics.median(ds[lonvar].values[pidx][idx])

                    surface_data['coords']['time']['data'] = np.append(surface_data['coords']['time']['data'], tm)
                    surface_data['data_vars']['deployment']['data'] = np.append(
                        surface_data['data_vars']['deployment']['data'],
                        deployment)
                    surface_data['data_vars']['depth']['data'] = np.append(surface_data['data_vars']['depth']['data'],
                                                                           np.nanmedian(ds.depth.values[pidx][idx]))
                    surface_data['data_vars']['lat']['data'] = np.append(surface_data['data_vars']['lat']['data'], gl_lat)
                    surface_data['data_vars']['lon']['data'] = np.append(surface_data['data_vars']['lon']['data'], gl_lon)
                    surface_data['data_vars']['pH_surface']['data'] = np.append(
                        surface_data['data_vars']['pH_surface']['data'],
                        pH_surface)
                    surface_data['data_vars']['omega_surface']['data'] = np.append(
                        surface_data['data_vars']['omega_surface']['data'],
                        omega_surface)
                    surface_data['data_vars']['temperature_surface']['data'] = np.append(
                        surface_data['data_vars']['temperature_surface']['data'],
                        temperature_surface)

                # grab bottom data
                idx = np.where(depth_interp[pidx] > max_profile_depth - 1)[0]

                # compare the glider depth to the global bathymetry file
                gl_lat = statistics.median(ds[latvar].values[pidx][idx])
                gl_lon = statistics.median(ds[lonvar].values[pidx][idx])
                profile_coords = [gl_lon, gl_lat]
                lat_idx = abs(bathy.lat.values - profile_coords[1]).argmin()
                lon_idx = abs(bathy.lon.values - profile_coords[0]).argmin()
                station_water_depth = -bathy.elevation[lat_idx, lon_idx].values

                # if the glider sample is within +/- 20% of the water column, append data to bottom dataset
                depth_threshold = [station_water_depth * .8, station_water_depth * 1.2]
                if np.logical_and(max_profile_depth > depth_threshold[0], max_profile_depth < depth_threshold[1]):
                    if np.sum(~np.isnan(ds[phvar].values[pidx][idx])) > 0:
                        tm = pd.to_datetime(ds.time.values[pidx][idx][0]).timestamp()
                        pH_bottom = np.round(np.nanmedian(ds[phvar].values[pidx][idx]), 4)
                        omega_bottom = np.round(np.nanmedian(arag[pidx][idx]), 2)
                        temperature_bottom = np.round(np.nanmedian(ds['temperature'].values[pidx][idx]), 2)
                        bottom_data['coords']['time']['data'] = np.append(bottom_data['coords']['time']['data'], tm)
                        bottom_data['data_vars']['deployment']['data'] = np.append(bottom_data['data_vars']['deployment']['data'],
                                                                            deployment)
                        bottom_data['data_vars']['depth']['data'] = np.append(bottom_data['data_vars']['depth']['data'],
                                                                       np.nanmedian(ds.depth.values[pidx][idx]))
                        bottom_data['data_vars']['bottom_depth']['data'] = np.append(bottom_data['data_vars']['bottom_depth']['data'],
                                                                              station_water_depth)
                        bottom_data['data_vars']['lat']['data'] = np.append(bottom_data['data_vars']['lat']['data'], gl_lat)
                        bottom_data['data_vars']['lon']['data'] = np.append(bottom_data['data_vars']['lon']['data'], gl_lon)
                        bottom_data['data_vars']['pH_bottom']['data'] = np.append(bottom_data['data_vars']['pH_bottom']['data'],
                                                                           pH_bottom)
                        bottom_data['data_vars']['omega_bottom']['data'] = np.append(bottom_data['data_vars']['omega_bottom']['data'],
                                                                              omega_bottom)
                        bottom_data['data_vars']['temperature_bottom']['data'] = np.append(
                            bottom_data['data_vars']['temperature_bottom']['data'],
                            temperature_bottom)

    # save bottom data as netcdf
    bottomds = xr.Dataset.from_dict(bottom_data)

    # add created time to global attrs
    datetime_format = '%Y-%m-%dT%H:%M:%SZ'
    created = dt.datetime.now(dt.UTC).strftime(datetime_format)  # creation time Timestamp
    time_start = dt.datetime.fromtimestamp(np.nanmin(bottomds.time.values), dt.UTC).strftime('%Y-%m-%d')
    time_end = dt.datetime.fromtimestamp(np.nanmax(bottomds.time.values), dt.UTC).strftime('%Y-%m-%d')
    start_yr = dt.datetime.fromtimestamp(np.nanmin(bottomds.time.values), dt.UTC).strftime('%Y')
    end_yr = dt.datetime.fromtimestamp(np.nanmax(bottomds.time.values), dt.UTC).strftime('%Y')

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

    global_attributes.update(bottomds.attrs)

    bottomds = bottomds.assign_attrs(global_attributes)
    bottomds = bottomds.sortby(bottomds.time)

    # Add compression to all variables
    encoding = {}
    for k in bottomds.data_vars:
        if k not in ['deployment']:
            encoding[k] = {'zlib': True, 'complevel': 1}

    encoding['time'] = dict(zlib=False, _FillValue=False, dtype=np.double)

    save_file = os.path.join(savedir, f'glider_based_bottom_OA_data_{start_yr}_{end_yr}.nc')
    bottomds.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims='time')

    # save surface data to netcdf
    surfaceds = xr.Dataset.from_dict(surface_data)

    # add created time to global attrs
    datetime_format = '%Y-%m-%dT%H:%M:%SZ'
    created = dt.datetime.now(dt.UTC).strftime(datetime_format)  # creation time Timestamp
    time_start = dt.datetime.fromtimestamp(np.nanmin(surfaceds.time.values), dt.UTC).strftime('%Y-%m-%d')
    time_end = dt.datetime.fromtimestamp(np.nanmax(surfaceds.time.values), dt.UTC).strftime('%Y-%m-%d')
    start_yr = dt.datetime.fromtimestamp(np.nanmin(surfaceds.time.values), dt.UTC).strftime('%Y')
    end_yr = dt.datetime.fromtimestamp(np.nanmax(surfaceds.time.values), dt.UTC).strftime('%Y')

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

    global_attributes.update(surfaceds.attrs)

    surfaceds = surfaceds.assign_attrs(global_attributes)
    surfaceds = surfaceds.sortby(surfaceds.time)

    # Add compression to all variables
    encoding = {}
    for k in surfaceds.data_vars:
        if k not in ['deployment']:
            encoding[k] = {'zlib': True, 'complevel': 1}

    encoding['time'] = dict(zlib=False, _FillValue=False, dtype=np.double)

    save_file = os.path.join(savedir, f'glider_based_surface_OA_data_{start_yr}_{end_yr}.nc')
    surfaceds.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims='time')


if __name__ == '__main__':
    rudir = '/Users/garzio/Documents/rucool/Saba/gliderdata/ru_ph_ncei_datasets'
    rudir2 = '/Users/garzio/Documents/rucool/Saba/gliderdata/ru_ph_ncei_datasets_prelim'
    sbudir = '/Users/garzio/Documents/rucool/Saba/gliderdata/sbu/from_DAC'
    umdir = '/Users/garzio/Documents/rucool/Saba/gliderdata/umaine/from_DAC'
    save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/output_nc'
    main(rudir, rudir2, sbudir, umdir, save_directory)
