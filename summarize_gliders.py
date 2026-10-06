#!/usr/bin/env python

"""
Author: Lori Garzio on 12/4/2024
Last modified: 10/6/2026
Export a .csv summary of glider deployments and dates
"""

import os
import numpy as np
import pandas as pd
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console

output_dir = os.path.join(os.getenv('HOME'), 'rucool/Saba/NOAA_SOE/data/output_nc')

# grab the dataset from ERDDAP
ru_server = 'https://rucool-sampling.marine.rutgers.edu/erddap'
dataset_id = 'vessel_glider_surface_bottom_OA_data_NES'
ds = cf.return_erddap_nc(ru_server, dataset_id)
ds = ds.swap_dims({'row': 'time'})

# subset for glider data
ds = ds.where(ds['collection_method'] == 'glider', drop=True)

deployments = np.unique(ds.cruise_deployment)
start_list = []
end_list = []
year_list = []
for d in deployments:
    idx = np.where(ds.cruise_deployment == str(d))[0]
    t0 = np.nanmin(ds.time[idx])
    tf = np.nanmax(ds.time[idx])
    start_list.append(pd.to_datetime(t0).strftime('%Y-%m-%d'))
    end_list.append(pd.to_datetime(tf).strftime('%Y-%m-%d'))
    year_list.append(pd.to_datetime(t0).strftime('%Y'))
    year_list.append(pd.to_datetime(tf).strftime('%Y'))

data = dict(
    glider_deployment=deployments,
    start=start_list,
    end=end_list
)

df = pd.DataFrame(data)
df.sort_values(by='start', inplace=True)

save_filename = f'glider_deployment_summary_{min(year_list)}_{max(year_list)}.csv'
save_file = os.path.join(output_dir, save_filename)
df.to_csv(save_file, index=False)