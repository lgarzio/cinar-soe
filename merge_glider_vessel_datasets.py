#!/usr/bin/env python

"""
Author: Lori Garzio on 8/5/2026
Last modified: 8/5/2026
Coming the glider- and vessel-based datasets together to create a single dataset per year
of surface and bottom pH and aragonite saturation state data for the U.S. Northeast Shelf.
"""

import os
import yaml
import glob
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
from collections import OrderedDict
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console
np.set_printoptions(suppress=True)



def main(glider_dir, vessel_file, savedir):
    home_dir = os.getenv('HOME')
    savedir = os.path.join(home_dir, savedir)
    os.makedirs(savedir, exist_ok=True)

    # get glider data
    glider_files = sorted(glob.glob(os.path.join(home_dir, glider_dir, '*.nc')))

    # combine all glider datasets
    datasets = [xr.open_dataset(f) for f in glider_files]
    gcombined = xr.concat(datasets, dim="time", data_vars="all", coords="minimal",
                          compat="override", join="outer", combine_attrs="override")
    gcombined = gcombined.sortby("time")

    vessel_ds = xr.open_dataset(os.path.join(home_dir, vessel_file))
    vessel_ds = vessel_ds.drop_vars(['data_source', 'obs_type', 'accession'], errors='ignore')

    # combine the glider and vessel datasets
    ds = xr.concat(
        [gcombined, vessel_ds],
        dim='time',
        data_vars='all',
        coords='minimal',
        compat='override',
        join='outer',
        combine_attrs='override',
    ).sortby('time')

    ds['obs'] = ('time', np.arange(ds.sizes['time']))
    ds['obs'].attrs['long_name'] = 'Observation number'

    # global attributes
    gafile = os.path.join(os.path.dirname(__file__), 'configs', 'global_attrs.yml')
    with open(gafile, "r") as file:
        global_attributes = yaml.safe_load(file)
    created = dt.datetime.now(dt.UTC).strftime('%Y-%m-%dT%H:%M')

    global_attributes['date_created'] = created
    global_attributes['date_modified'] = created

    # encoding for variables
    encoding = {}
    for k in ds.data_vars:
        if k not in ['cruise_deployment']:
            encoding[k] = {'zlib': True, 'complevel': 1}

    encoding['time'] = dict(zlib=False, _FillValue=None, dtype=np.double)

    # separate the datasets into yearly files
    for year in np.unique(ds['time.year']):
        ds_year = ds.sel(time=ds['time.year'] == year)

        # add time coverage start and end to global attributes
        time_start = pd.to_datetime(ds_year.time.values.min()).strftime('%Y-%m-%d')
        time_end = pd.to_datetime(ds_year.time.values.max()).strftime('%Y-%m-%d')
        global_attributes['time_coverage_start'] = time_start
        global_attributes['time_coverage_end'] = time_end

        ds_year = ds_year.assign_attrs(global_attributes)
        ds_year = ds_year.sortby(ds_year.time)

        ds_year['time'] = xr.DataArray(
            (pd.to_datetime(ds_year.time.values) - pd.Timestamp('1970-01-01')).astype('timedelta64[s]').astype(np.int64),
            dims=('time',),
            coords={'time': ds_year.time.values},
            attrs={'units': 'seconds since 1970-01-01 00:00:00', 'standard_name': 'time', 'long_name': 'time'},
        )

        ds_year = ds_year.swap_dims({'time': 'obs'})
        save_file = os.path.join(savedir, f'{str(year)}_surface_bottom_OA_data.nc')
        ds_year.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims='time')


if __name__ == '__main__':
    glider_files = 'rucool/Saba/NOAA_SOE/data/output_nc/gliders'
    vessel_file = 'rucool/Saba/NOAA_SOE/data/output_nc/vessel_based_OA_data_2004_2024.nc'
    save_directory = 'rucool/Saba/NOAA_SOE/data/output_nc/merged'
    main(glider_files, vessel_file, save_directory)
