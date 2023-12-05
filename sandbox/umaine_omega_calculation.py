import os
import numpy as np
import pandas as pd
import xarray as xr
import functions.common as cf

# calculate TA and aragonite saturation state using equations IV and VI from Table 3 of McGarry et al 2020

f = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/nc_files/um_242-20210630T1916-profile-sci-delayed_shifted.nc'
ds = xr.open_dataset(f)

temp = ds.temperature.values
sal = ds.salinity.values
oxy = ds.oxygen_concentration.values

# calculate TA and add to dataset
# normalize temperature and salinity to the McGarry et al 2020 data from Table 3
tempn = (temp - 13.20) / 5.92
saln = (sal - 34.40) / 1.49
ta = 2289 + (0.758 * tempn) + (69.2 * saln)


attrs = {
        'ancillary_variables': 'salinity',
        'observation_type': 'calculated',
        'units': 'umol/kg',
        'long_name': 'Total Alkalinity',
        'comment': 'Calculated from temperature and salinity using equation IV from Table 3 in McGarry et al 2020:'
                   'https://doi.org/10.1029/2020JC016480'
    }
da = xr.DataArray(ta, coords=ds['salinity'].coords, dims=ds['salinity'].dims,
                  name='total_alkalinity', attrs=attrs)
ds['total_alkalinity'] = da

# # calculate omega and add to dataset
# omega = 0.02065 + (0.00635 * temp)
# omega2 = (0.0681 * sal) + (0.000101 * temp * sal)
#
# attrs = {
#         'ancillary_variables': 'temperature salinity',
#         'observation_type': 'calculated',
#         'units': '1',
#         'long_name': 'Aragonite Saturation State',
#         'comment': 'Calculated from salinity and temperature using equation VI from Table 3 in McGarry et al 2020:'
#                    'https://doi.org/10.1029/2020JC016480. (0.0681 * sal) + (0.000101 * temp * sal)'
#     }
# da = xr.DataArray(omega2, coords=ds['salinity'].coords, dims=ds['salinity'].dims,
#                   name='saturation_aragonite', attrs=attrs)
# ds['saturation_aragonite'] = da

omega_arag, pco2, revelle = cf.run_co2sys_ta_ph(ds.total_alkalinity.values,
                                                ds.ph_total_shifted.values,
                                                ds.salinity.values,
                                                ds.temperature.values,
                                                ds.pressure.values)

attrs = {
        'ancillary_variables': 'total_alkalinity ph_total_shifted temperature salinity pressure',
        'observation_type': 'calculated',
        'units': '1',
        'long_name': 'Aragonite Saturation State',
        'comment': 'Calculated using the PyCO2SYS function with inputs of Total Alkalinity, pH, salinity, temperature, and pressure'
    }
da = xr.DataArray(omega_arag, coords=ds['salinity'].coords, dims=ds['salinity'].dims,
                  name='saturation_aragonite', attrs=attrs)
ds['saturation_aragonite'] = da

fname = '{}_final.nc'.format(f.split('.')[0])
ds.to_netcdf(fname)