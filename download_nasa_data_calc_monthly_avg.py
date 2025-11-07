#!/usr/bin/env python

"""
Author: Lori Garzio on 11/15/2024
Last modified: 11/7/2025
Download NASA 4km resolution SST data from opendap: monthly and daily files.
https://oceandata.sci.gsfc.nasa.gov/l3/ data viewer and the backend opendap server for the SNPP VIIRS L3 dataset:
is https://oceandata.sci.gsfc.nasa.gov/opendap/hyrax/VIIRS/L3SMI/contents.html 
example monthly file: http://oceandata.sci.gsfc.nasa.gov/opendap/VIIRS/L3SMI/2012/0301/SNPP_VIIRS.20120301_20120331.L3m.MO.SST3.sst_triple.4km.nc
example daily file: http://oceandata.sci.gsfc.nasa.gov/opendap/VIIRS/L3SMI/2024/0604/SNPP_VIIRS.20240604.L3m.DAY.SST3.sst_triple.4km.nc
Monthly SST files were downloaded from March 2012 - May 2024
Monthly files are no longer available so download daily files for June 2024+ then calculate monthly averages.

Also download monthly NASA chla data
example monthly file: http://oceandata.sci.gsfc.nasa.gov/opendap/VIIRS/L3SMI/2012/0301/SNPP_VIIRS.20120301_20120331.L3m.MO.CHL.chlor_a.4km.nc
"""

import os
import glob
import xarray as xr
import numpy as np
import pandas as pd
import calendar
import datetime as dt


def main(extent, frequency, start, end, v, savedir):
    savedir = os.path.join(savedir, v.split('.')[-1])
    os.makedirs(savedir, exist_ok=True)
    if frequency == 'monthly':
        savedirmonthly = os.path.join(savedir, 'opendap_monthly')
        os.makedirs(savedirmonthly, exist_ok=True)

        # download monthly files
        months = pd.date_range(start, end, freq='MS')
        for sday in months:
            year = sday.year
            month = sday.month
            mm = str(month).zfill(2)
            start_ymd = sday.strftime('%Y%m%d')
            end_day = calendar.monthrange(year, month)[1]
            end_ymd = (sday + dt.timedelta(days=(end_day - 1))).strftime('%Y%m%d')

            dirname = f'http://oceandata.sci.gsfc.nasa.gov/opendap/VIIRS/L3SMI/{year}/{mm}01'
            fname = f'SNPP_VIIRS.{start_ymd}_{end_ymd}.L3m.MO.{v}.4km.nc'

            try:
                ds = xr.open_dataset(os.path.join(dirname, fname))
                lati = np.where(np.logical_and(ds.lat.values > extent[2], ds.lat.values < extent[3]))[0]
                loni = np.where(np.logical_and(ds.lon.values > extent[0], ds.lon.values < extent[1]))[0]
                ds = ds.isel(lat=lati, lon=loni)

                ds.to_netcdf(os.path.join(savedirmonthly, fname))
            except OSError:
                print(f'file not found {os.path.join(dirname, fname)}')
                continue

    elif frequency == 'daily':
        # download daily files
        savedirdaily = os.path.join(savedir, 'opendap_daily')

        days = pd.date_range(start, end, freq='D')
        for day in days:
            year = day.year
            month = day.month
            mm = str(month).zfill(2)
            savedirmonth = os.path.join(savedirdaily, f'{year}{mm}')
            os.makedirs(savedirmonth, exist_ok=True)
            daynum = day.day
            dd = str(daynum).zfill(2)
            start_ymd = day.strftime('%Y%m%d')

            #f = 'http://oceandata.sci.gsfc.nasa.gov/opendap/VIIRS/L3SMI/2024/0604/SNPP_VIIRS.20240604.L3m.DAY.SST3.sst_triple.4km.nc'
            dirname = f'http://oceandata.sci.gsfc.nasa.gov/opendap/VIIRS/L3SMI/{year}/{mm}{dd}'
            fname = f'SNPP_VIIRS.{start_ymd}.L3m.DAY.SST3.sst_triple.4km.nc'

            try:
                ds = xr.open_dataset(os.path.join(dirname, fname))
                lati = np.where(np.logical_and(ds.lat.values > extent[2], ds.lat.values < extent[3]))[0]
                loni = np.where(np.logical_and(ds.lon.values > extent[0], ds.lon.values < extent[1]))[0]
                ds = ds.isel(lat=lati, lon=loni)

                ds.to_netcdf(os.path.join(savedirmonth, fname))
            except OSError:
                print(f'file not found {os.path.join(dirname, fname)}')
                continue

        # Calculate monthly averages from the daily files
        months = pd.date_range(start, end, freq='MS')
        for sday in months:
            year = sday.year
            month = sday.month
            mm = str(month).zfill(2)

            filedir = os.path.join(savedir, 'opendap_daily', f'{year}{mm}')
            sfile = os.path.join(savedir, 'opendap_monthly', f'SNPP_VIIRS.{year}{mm}01_{year}{mm}{calendar.monthrange(year, month)[1]}.L3m.MO.{v}.4km.nc')
            file_list = glob.glob(os.path.join(filedir, '*.nc'))
            for i, f in enumerate(file_list):
                ds = xr.open_dataset(f)
                ds = ds.drop_vars(['palette', 'qual_sst_triple'])
                if i == 0:
                    data = [ds[v.split('.')[-1]].values]
                else:
                    data = np.concatenate([data, [ds[v.split('.')[-1]].values]])

            data_avg = np.nanmean(data, axis=0)

            ds[v.split('.')[-1]].values = data_avg
            ds.attrs['comment'] = f'Monthy average: {year}-{mm}'
            ds[v.split('.')[-1]].attrs['comment'] = f'Monthy average: {year}-{mm}'

            ds.to_netcdf(sfile)


if __name__ == '__main__':
    extent = [-79, -64, 34, 46]
    frequency = 'monthly'  # 'monthly' or 'daily' - download daily files then calculate monthly averages
    startdate = '20120201'
    enddate = '20251031'
    variable = 'CHL.chlor_a'  # 'SST3.sst_triple' or 'CHL.chlor_a'
    # save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/nasa_data'
    save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE/data/nasa_data'
    main(extent, frequency, startdate, enddate, variable, save_directory)