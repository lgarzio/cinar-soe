#!/usr/bin/env python

"""
Author: Lori Garzio on 8/5/2026
Last modified: 10/6/2026
Combine the glider- and vessel-based datasets together to create a single dataset per year
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

    # attributes config file
    attrsfile = os.path.join(os.path.dirname(__file__), 'configs', 'merged_attrs.yml')
    with open(attrsfile, "r") as file:
        attributes_config = yaml.safe_load(file)
    created = dt.datetime.now(dt.UTC).strftime('%Y-%m-%dT%H:%M')
    attributes_config['global_attrs']['date_created'] = created
    attributes_config['global_attrs']['date_modified'] = created
    
    # get glider data
    glider_files = sorted(glob.glob(os.path.join(home_dir, glider_dir, '*.nc')))

    # combine all glider datasets
    datasets = [xr.open_dataset(f) for f in glider_files]
    gcombined = xr.concat(datasets, dim="time", data_vars="all", coords="minimal",
                          compat="override", join="outer", combine_attrs="override")
    gcombined = gcombined.sortby("time")
    
    gcombined["collection_method"] = (
        "time",
        np.full(gcombined.sizes["time"], "glider", dtype="<U6"),
    )

    vessel_ds = xr.open_dataset(os.path.join(home_dir, vessel_file))
    vessel_ds = vessel_ds.drop_vars(['data_source', 'obs_type', 'accession'], errors='ignore')

    vessel_ds["collection_method"] = (
            "time",
            np.full(vessel_ds.sizes["time"], "vessel", dtype="<U6"),
        )

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

    # add variable attributes from the config file
    for var in ds.data_vars:
        if var in attributes_config['variable_attrs']:
            ds[var] = ds[var].assign_attrs(attributes_config['variable_attrs'][var])

    # encoding for variables
    encoding = {}
    for k in ds.data_vars:
        if k not in ['cruise_deployment', 'collection_method']:
            encoding[k] = {'zlib': True, 'complevel': 1}

    encoding['time'] = dict(zlib=False, _FillValue=None, dtype=np.double)

    # separate the datasets into yearly files
    for year in np.unique(ds['time.year']):
        ds_year = ds.sel(time=ds['time.year'] == year)

        # add time coverage start and end to global attributes
        time_start = pd.to_datetime(ds_year.time.values.min()).strftime('%Y-%m-%d')
        time_end = pd.to_datetime(ds_year.time.values.max()).strftime('%Y-%m-%d')
        attributes_config['global_attrs']['time_coverage_start'] = time_start
        attributes_config['global_attrs']['time_coverage_end'] = time_end

        ds_year = ds_year.assign_attrs(attributes_config['global_attrs'])
        ds_year = ds_year.sortby(ds_year.time)

        ds_year['time'] = xr.DataArray(
            (pd.to_datetime(ds_year.time.values) - pd.Timestamp('1970-01-01')).astype('timedelta64[s]').astype(np.int64),
            dims=('time',),
            coords={'time': ds_year.time.values},
            attrs={'units': 'seconds since 1970-01-01 00:00:00', 'standard_name': 'time', 'long_name': 'time'},
        )

        #ds_year = ds_year.swap_dims({'time': 'obs'})
        save_file = os.path.join(savedir, f'glider_vessel_surface_bottom_OA_data_{str(year)}.nc')
        ds_year.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims='time')


if __name__ == '__main__':
    glider_files = 'rucool/Saba/NOAA_SOE/data/output_nc/gliders'
    vessel_file = 'rucool/Saba/NOAA_SOE/data/output_nc/vessel_based_OA_data_2004_2024.nc'
    save_directory = 'rucool/Saba/NOAA_SOE/data/output_nc/merged'
    main(glider_files, vessel_file, save_directory)
