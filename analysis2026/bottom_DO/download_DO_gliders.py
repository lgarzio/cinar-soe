#!/usr/bin/env python

"""
Author: Lori Garzio on 10/22/2025
Last modified: 10/22/2025
Download glider datasets from the IOOS Glider DAC
https://gliders.ioos.us/erddap/index.html in the MAB that contain
dissolved oxygen for a user-specified time period
"""

import os
import glob
import datetime as dt
import numpy as np
import pandas as pd
import xarray as xr
from erddapy import ERDDAP
from collections import Counter
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console
np.set_printoptions(suppress=True)


def return_dataset_variables(server, dataset_id):
    e = ERDDAP(server=server,
               protocol='tabledap',
               response='nc')
    var_dict = e._get_variables(dataset_id=dataset_id)

    return list(var_dict.keys())


def return_erddap_nc(server, dataset_id, variables=None, constraints=None):
    e = ERDDAP(server=server,
               protocol='tabledap',
               response='nc')
                      
    e.dataset_id = dataset_id
    if constraints:
        e.constraints = constraints
    if variables:
        e.variables = variables
    
    ds = e.to_xarray(requests_kwargs={"timeout": 600})  # increase timeout to 10 minutes
    ds = ds.sortby(ds.time)
    return ds


def return_glider_ids(server, kwargs):
    """
    Searches an ERDDAP server for datasets and returns dataset IDs
    :param server: e.g. 'https://gliders.ioos.us/erddap'
    :param kwargs: dictionary containing coordinate and time limits
    :return: array containing dataset IDs
    """
    e = ERDDAP(server=server)
    search_url = e.get_search_url(response='csv', **kwargs)
    search = pd.read_csv(search_url)
    ds_ids = search['Dataset ID'].values

    return ds_ids


def main(savedir, min_time, max_time):
    # grab all gliders in the MAB region for 2025 from the IOOS Glider DAC
    kwargs = dict()
    kwargs['min_time'] = min_time
    kwargs['max_time'] = max_time
    kwargs['min_lon'] = -78
    kwargs['max_lon'] = -65
    kwargs['min_lat'] = 35
    kwargs['max_lat'] = 45

    #server = 'https://gliders.ioos.us/erddap'
    server = 'https://slocum-data.marine.rutgers.edu/erddap'
    glider_ids = return_glider_ids(server, kwargs)

    # for each glider deployment, figure out if dissolved oxygen data are available
    glider_do_ids = dict()
    for gi in glider_ids:
        if '-trajectory-raw-' in gi:
            continue
        ds_vars = return_dataset_variables(server, gi)
        oxygen_vars = [var for var in ds_vars if 'oxygen' in var.lower()]
        #print(oxygen_vars)
        if len(oxygen_vars) > 0:
            glider_do_ids[gi] = oxygen_vars

    # if there are rt and delayed mode datasets, grab only the delayed mode
    check_ids = [gi.replace('-delayed', '') for gi in list(glider_do_ids.keys())]
    duplicates = [item for item, count in Counter(check_ids).items() if count > 1]
    glider_do_ids_final = {key: value for key, value in glider_do_ids.items() if key not in duplicates}

    # download the datasets to your local machine
    #vars = ['time', 'latitude', 'longitude', 'depth', 'profile_id', 'temperature']
    vars = ['time', 'latitude', 'longitude', 'depth', 'profile_time', 'profile_lon', 'profile_lat', 'temperature']
    skip = ['maracoos_04-20241203T1457-profile-sci-delayed', 'maracoos_05-20250404T1319-profile-sci-delayed',
            'maracoos_06-20250116T1522-profile-sci-rt', 'maracoos_06-20250520T1554-profile-sci-rt',
            'maracoos_06-20250929T1630-profile-sci-rt', 'ru32-20250716T1541-profile-sci-delayed']
    for gi, oxygen_vars in glider_do_ids_final.items():
        if gi in skip:
            print(f'skipping {gi}')
            continue
        print(f'downloading {gi}')
        
        downloadvars = vars + oxygen_vars
        
        kwargs = dict()
        kwargs['variables'] = downloadvars
        ds = return_erddap_nc(server, gi, **kwargs)

        fname = f'{gi}.nc'
        ds.to_netcdf(os.path.join(savedir, fname))


if __name__ == '__main__':
    save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/dissolved_oxygen_gliders/files/rutgers'
    start = "2025-01-01T00:00:00Z"
    end = "2025-10-22T23:59:59Z"
    main(save_directory, start, end)
