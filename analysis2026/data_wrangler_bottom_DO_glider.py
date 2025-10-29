#!/usr/bin/env python

"""
Author: Lori Garzio on 10/23/2025
Last modified: 10/29/2025
Grab bottom-water DO from glider datasets downloaded from the IOOS Glider DAC 
https://gliders.ioos.us/erddap/index.html (Stony Brook gliders)
or the Rutgers ERDDAP server https://slocum-data.marine.rutgers.edu/erddap/index.html
using download_DO_gliders.py and export as NetCDF.
Files were downloaded on 10/22/2025
"""

import os
import glob
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
from collections import OrderedDict
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console
np.set_printoptions(suppress=True)


def interpolate_var(da):
    df = da.to_dataframe()

    # drop duplicated columns
    df = df.loc[:, ~df.columns.duplicated()]

    var_interp = df[da.name].interpolate(method='linear', limit_direction='both', limit=2).values
    return var_interp


def main(fdir, lon_bounds, lat_bounds, yr):
    extent = [lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]]
    # duplicate the lon/lat bounds to make a closed polygon for shapely
    lon_bounds.append(lon_bounds[1])
    lon_bounds.append(lon_bounds[0])
    lat_bounds = [lat_bounds[0], lat_bounds[0], lat_bounds[1],lat_bounds[1]]

    bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))

    # initialize dictionary to append bottom DO data from glider deployments
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
            "comment": f"Synthesis of bottom DO data from glider-based measurements downloaded from the IOOS Glider DAC"
                       "https://gliders.ioos.us/erddap/index.html for {yr} that were spatially limited to the U.S. Northeast Shelf."
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
            "DO_bottom": {
                "dims": "time",
                "data": np.array([], dtype='float32'),
                "attrs": {
                    "units": "mg L-1",
                    "long_name": "Bottom Dissolved Oxygen",
                    "description": "Bottom DO is the the median of the values within the deepest 1m of a glider "
                                   "profile, provided the sampling depth was no shallower than the bottom 20% of total "
                                   "water column depth."
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

    oxygen_varnames = ['dissolved_oxygen', 'oxy4_oxygen', 'oxygen', 'oxy3835_wphase_oxygen', 'oxygen_concentration_shifted_mgL']

    for gf in sorted(glob.glob(os.path.join(fdir, '*.nc'))):
        print(gf)
        ds = xr.open_dataset(gf)
        ds = ds.drop_dims('trajectory')
        deployment = ds.title.split(' ')[0]
        depth_interp = interpolate_var(ds.depth)

        # figure out what the dissolved oxygen variable is called in this dataset
        oxyvar = set(oxygen_varnames).intersection(set(ds.data_vars))

        if len(oxyvar) == 0:
            raise ValueError(f'No recognized dissolved oxygen variable found in {deployment}')
        elif len(oxyvar) > 1:
            raise ValueError(f'Multiple recognized dissolved oxygen variables found in {deployment}')
        else:
            oxyvar = list(oxyvar)[0]

        #print(f'{oxyvar} ({ds[oxyvar].units})')
        # convert DO to mg/L
        if ds[oxyvar].units in['umol kg-1', 'micromol kg-1', 'micromoles L-1', 'umol kg-3']:  # umol kg-3 is a typo it should be umol L-1
            oxydata = ds[oxyvar].values * 32 / 1000  # convert from umol/L aka roughly umol/kg to mg/L
        elif ds[oxyvar].units in ['mg L-1']:
            oxydata = ds[oxyvar].values
        else:
            raise ValueError(f'Unrecognized DO units ({ds[oxyvar].units}) in {deployment}')

        # for each glider, grab the data in the bottom 1m of each profile when the glider is within +/- 20% of the water
        # column depth (determined by comparison to global bathymetry data)
        start_idx = 0
        for pid in ds.profile.values:
            pds = ds.sel(profile=pid)

            # figure out which index corresponds to which profile based on rowSize
            pidx = np.arange(start_idx, start_idx + pds.rowSize.values)
            start_idx += pds.rowSize.values
            
            profile_year = pd.to_datetime(pds.time.values).year
            gllat = pds.latitude.values
            gllon = pds.longitude.values
            
            if isinstance(profile_year, pd.core.indexes.base.Index):
                profile_year = pd.to_datetime(ds.profile_time.values[pidx][0]).year
                gllat = ds.profile_lat.values[pidx][0]
                gllon = ds.profile_lon.values[pidx][0]
            
            # make sure the profile was collected during the specified year
            if profile_year != yr:
                continue
            
            # make sure the profile is within the defined region
            if Polygon(list(zip(lon_bounds, lat_bounds))).contains(Point(gllon, gllat)):
                max_profile_depth = np.nanmax(ds.depth.values[pidx])
                if max_profile_depth > 10:

                    # grab bottom data
                    idx = np.where(depth_interp[pidx] > max_profile_depth - 1)[0]

                    # compare the glider depth to the global bathymetry file
                    profile_coords = [gllon, gllat]
                    lat_idx = abs(bathy.lat.values - profile_coords[1]).argmin()
                    lon_idx = abs(bathy.lon.values - profile_coords[0]).argmin()
                    station_water_depth = -bathy.elevation[lat_idx, lon_idx].values

                    # if the glider sample is within +/- 20% of the water column, append data to bottom dataset
                    depth_threshold = [station_water_depth * .8, station_water_depth * 1.2]
                    if np.logical_and(max_profile_depth > depth_threshold[0], max_profile_depth < depth_threshold[1]):
                        if np.sum(~np.isnan(oxydata[pidx][idx])) > 0:
                            try:
                                tm = pd.to_datetime(pds.time.values).timestamp()
                            except AttributeError:
                                tm = pd.to_datetime(ds.profile_time.values[pidx][0]).timestamp()

                            DO_bottom = np.round(np.nanmedian(oxydata[pidx][idx]), 4)
                            if DO_bottom == 0.0:
                                continue
                            if deployment == 'ru33-20250903T1642':  # get rid of one outlier
                                if DO_bottom < 3:
                                    continue

                            temperature_bottom = np.round(np.nanmedian(ds['temperature'].values[pidx][idx]), 2)
                            bottom_data['coords']['time']['data'] = np.append(bottom_data['coords']['time']['data'], tm)
                            bottom_data['data_vars']['deployment']['data'] = np.append(bottom_data['data_vars']['deployment']['data'],
                                                                                deployment)
                            bottom_data['data_vars']['depth']['data'] = np.append(bottom_data['data_vars']['depth']['data'],
                                                                        np.nanmedian(depth_interp[pidx][idx]))
                            bottom_data['data_vars']['bottom_depth']['data'] = np.append(bottom_data['data_vars']['bottom_depth']['data'],
                                                                                station_water_depth)
                            bottom_data['data_vars']['lat']['data'] = np.append(bottom_data['data_vars']['lat']['data'], gllat)
                            bottom_data['data_vars']['lon']['data'] = np.append(bottom_data['data_vars']['lon']['data'], gllon)
                            bottom_data['data_vars']['DO_bottom']['data'] = np.append(bottom_data['data_vars']['DO_bottom']['data'],
                                                                            DO_bottom)
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
    print(start_yr, end_yr)

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

    save_file = os.path.join(os.path.dirname(fdir), f'glider_based_bottom_DO_data_{yr}.nc')
    bottomds.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims='time')


if __name__ == '__main__':
    file_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/dissolved_oxygen_gliders/files'
    lon_extent = [-76, -68]  # longitude boundaries for grabbing data in the MAB
    lat_extent = [38, 42]  # latitude boundaries for grabbing data in the MAB
    year = 2025
    main(file_directory, lon_extent, lat_extent, year)
