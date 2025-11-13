#!/usr/bin/env python

"""
Author: Lori Garzio on 12/4/2024
Last modified: 11/13/2025
Export a .csv summary of glider deployments and dates
"""

import numpy as np
import pandas as pd
import xarray as xr
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console


file = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/output_nc/glider_based_surface_OA_data_2019_2025.nc'
save_file = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/output_nc/glider_deployment_summary_2019_2025.csv'

ds = xr.open_dataset(file)

deployments = np.unique(ds.deployment)
start_list = []
end_list = []
for d in deployments:
    idx = np.where(ds.deployment == str(d))[0]
    t0 = np.nanmin(ds.time[idx])
    tf = np.nanmax(ds.time[idx])
    start_list.append(pd.to_datetime(t0).strftime('%Y-%m-%d'))
    end_list.append(pd.to_datetime(tf).strftime('%Y-%m-%d'))

data = dict(
    glider_deployment=deployments,
    start=start_list,
    end=end_list
)

df = pd.DataFrame(data)
df.sort_values(by='start', inplace=True)
df.to_csv(save_file, index=False)