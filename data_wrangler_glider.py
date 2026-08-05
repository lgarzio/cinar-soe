#!/usr/bin/env python

"""
Author: Lori Garzio on 10/23/2024
Last modified: 8/5/2026
Grab bottom- and surface-water pH and omega data from glider datasets and export as NetCDF.
Datasets are available on the IOOS Glider DAC ERDDAP server https://gliders.ioos.us/erddap/index.html
as well as the NCEI OCADS data portal https://www.ncei.noaa.gov/products/ocean-carbon-acidification-data-system.
"""

import os
import yaml
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


def interpolate_var(ds, varname):
    df = ds.to_dataframe()
    var_interp = df[varname].interpolate(method='linear', limit_direction='both', limit=2).values
    newvarname = f'{varname}_interpolated'
    ds[newvarname] = xr.DataArray(var_interp, coords=ds[varname].coords, dims=ds[varname].dims, name=newvarname)


def main(filedir, savedir):
    home_dir = os.getenv('HOME')
    savedir = os.path.join(home_dir, savedir)
    bathymetry = os.path.join(home_dir, 'rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc')
    extent = [-78, -65, 35, 45]
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))

    # get variable configuration file
    cfile = os.path.join(os.path.dirname(__file__), 'configs', 'glider_vars.yml')

    # get glider data
    glider_files = sorted(glob.glob(os.path.join(home_dir, filedir, '*.nc')))

    for gf in glider_files:
        print(gf)

        # initialize dictionary to append surface and bottom pH and aragonite data from glider deployment
        data = {
            "attrs": {
                "comment": "Synthesis of surface and bottom pH and aragonite saturation state data from glider-based measurements "
                            "that were spatially limited to the U.S. Northeast Shelf.",
            },
            "dims": "time",
            "coords": {},
            "data_vars": {}
        }
        
        # add variables from the config file to the data dictionary
        with open(cfile, "r") as file:
            config_dict = yaml.safe_load(file)
        for varname, var_dict in config_dict.items():
            dtype = var_dict.pop('dtype')
            var_dict['data'] = np.array([], dtype=dtype)
            if varname == 'time':
                data['coords'][varname] = var_dict
            else:
                data['data_vars'][varname] = var_dict

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
            try:
                latvar = 'lat'
                lonvar = 'lon'
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

        # interpolate temperature
        interpolate_var(ds, 'temperature')

        # calculate TA and omega for sbu gliders
        if glider in ['sbu01', 'sbu02']:
            interpolate_var(ds, 'salinity')
            interpolate_var(ds, 'temperature')
            interpolate_var(ds, 'pressure')
            season = str(ds['time.season'].values[0])
            ta = cf.calc_ta_nyb(season, ds['salinity_interpolated'].values)

            # when pH is nan, TA is nan
            idx = np.isnan(ds[phvar].values)
            ta[idx] = np.nan

            # run CO2SYS to calculate omega
            arag, pco2, revelle = cf.run_co2sys_ta_ph(ta,
                                                      ds[phvar].values,
                                                      ds['salinity_interpolated'].values,
                                                      ds['temperature_interpolated'].values,
                                                      ds['pressure_interpolated'].values)
            ds['aragonite_saturation_state'] = xr.DataArray(arag, coords=ds[phvar].coords, dims=ds[phvar].dims, name='aragonite_saturation_state')

        # apply QC to umaine gliders for some issues I noticed
        if glider in ['um_242']:
            # TA
            idx = np.where(ds.total_alkalinity.values < 1900)[0]
            ds[phvar].values[idx] = np.nan
            ds['aragonite_saturation_state'].values[idx] = np.nan

            # omega
            idx = np.where(ds['aragonite_saturation_state'].values > 5)[0]
            ds[phvar].values[idx] = np.nan
            ds['aragonite_saturation_state'].values[idx] = np.nan

            # pH
            idx = np.where(ds[phvar] > 8.6)[0]
            ds[phvar].values[idx] = np.nan
            ds['aragonite_saturation_state'].values[idx] = np.nan

        try:
            ds.depth_interpolated.values
        except AttributeError:
            interpolate_var(ds, 'depth')

        df_vars = ['time','depth_interpolated', phvar, 'aragonite_saturation_state', 'temperature_interpolated', latvar, lonvar]

        # for each glider, grab the data in the bottom 1m of each profile when the glider is within +/- 20% of the water
        # column depth (determined by comparison to global bathymetry data)
        # also grab the surface data (2-4m depth)
        for pid in profileid:
            ds_sel = ds.where(ds[profile_identifier] == pid, drop=True)

            # check for any valid pH data in the profile, if not skip to the next profile
            if np.sum(~np.isnan(ds_sel[phvar].values)) == 0:
                continue

            max_profile_depth = np.nanmax(ds_sel.depth_interpolated.values)
            if max_profile_depth > 10:
                df_sel = ds_sel[df_vars].to_dataframe().reset_index()
                df_sel = df_sel.set_index('time')
                tm = pd.to_datetime(df_sel.index).mean().timestamp()
                gl_lat = np.nanmedian(df_sel[latvar].values)
                gl_lon = np.nanmedian(df_sel[lonvar].values)

                # grab suface data (2-4m depth)
                df_sel_surface = df_sel[df_sel['depth_interpolated'].between(2, 4, inclusive='both')]
                depth_surface = np.nanmedian(df_sel_surface['depth_interpolated'].values)
                pH_surface = np.nanmedian(df_sel_surface[phvar].values)
                omega_surface = np.nanmedian(df_sel_surface['aragonite_saturation_state'].values)
                temp_surface = np.nanmedian(df_sel_surface['temperature_interpolated'].values)

                # grab bottom data (max depth - 1m)
                # compare the glider depth to the global bathymetry file
                profile_coords = [gl_lon, gl_lat]
                lat_idx = abs(bathy.lat.values - profile_coords[1]).argmin()
                lon_idx = abs(bathy.lon.values - profile_coords[0]).argmin()
                station_water_depth = -bathy.elevation[lat_idx, lon_idx].values

                # if the glider sample is within 20% of the depth of the water column, append data to bottom dataset
                if max_profile_depth > station_water_depth * .8:
                    df_sel_bottom = df_sel[df_sel['depth_interpolated'] > max_profile_depth - 1]
                    depth_bottom = np.nanmedian(df_sel_bottom['depth_interpolated'].values)
                    pH_bottom = np.nanmedian(df_sel_bottom[phvar].values)
                    omega_bottom = np.nanmedian(df_sel_bottom['aragonite_saturation_state'].values)
                    temp_bottom = np.nanmedian(df_sel_bottom['temperature_interpolated'].values)
                else:
                    depth_bottom = np.nan
                    pH_bottom = np.nan
                    omega_bottom = np.nan
                    temp_bottom = np.nan

                if np.isnan(pH_surface) and np.isnan(pH_bottom):
                    print(f'No valid pH surface or bottom data for profile {pid} in deployment {deployment}, skipping...')
                    continue

                # add data to dictionary
                data['coords']['time']['data'] = np.append(data['coords']['time']['data'], tm)

                data['data_vars']['cruise_deployment']['data'] = np.append(data['data_vars']['cruise_deployment']['data'], deployment)
                
                data['data_vars']['depth_surface']['data'] = np.append(data['data_vars']['depth_surface']['data'], depth_surface)
                data['data_vars']['depth_bottom']['data'] = np.append(data['data_vars']['depth_bottom']['data'], depth_bottom)
                data['data_vars']['station_water_depth']['data'] = np.append(data['data_vars']['station_water_depth']['data'], station_water_depth)
                data['data_vars']['lat']['data'] = np.append(data['data_vars']['lat']['data'], gl_lat)
                data['data_vars']['lon']['data'] = np.append(data['data_vars']['lon']['data'], gl_lon)

                data['data_vars']['pH_surface']['data'] = np.append(data['data_vars']['pH_surface']['data'], pH_surface)
                data['data_vars']['omega_surface']['data'] = np.append(data['data_vars']['omega_surface']['data'], omega_surface)
                data['data_vars']['temperature_surface']['data'] = np.append(data['data_vars']['temperature_surface']['data'], temp_surface)
                data['data_vars']['pH_bottom']['data'] = np.append(data['data_vars']['pH_bottom']['data'], pH_bottom)
                data['data_vars']['omega_bottom']['data'] = np.append(data['data_vars']['omega_bottom']['data'], omega_bottom)
                data['data_vars']['temperature_bottom']['data'] = np.append(data['data_vars']['temperature_bottom']['data'], temp_bottom)

        # save data as netcdf
        ds = xr.Dataset.from_dict(data)

        # add created time to global attrs
        datetime_format = '%Y-%m-%dT%H:%M:%SZ'
        created = dt.datetime.now(dt.UTC).strftime(datetime_format)  # creation time Timestamp
        time_start = dt.datetime.fromtimestamp(np.nanmin(ds.time.values), dt.UTC).strftime('%Y-%m-%d')
        time_end = dt.datetime.fromtimestamp(np.nanmax(ds.time.values), dt.UTC).strftime('%Y-%m-%d')
        start_yr = dt.datetime.fromtimestamp(np.nanmin(ds.time.values), dt.UTC).strftime('%Y')
        end_yr = dt.datetime.fromtimestamp(np.nanmax(ds.time.values), dt.UTC).strftime('%Y')

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

        global_attributes.update(ds.attrs)

        ds = ds.assign_attrs(global_attributes)
        ds = ds.sortby(ds.time)

        # Add compression to all variables
        encoding = {}
        for k in ds.data_vars:
            if k not in ['cruise_deployment']:
                encoding[k] = {'zlib': True, 'complevel': 1}

        encoding['time'] = dict(zlib=False, _FillValue=False, dtype=np.double)

        save_file = os.path.join(savedir, f'{deployment}_surface_bottom_OA_data.nc')
        ds.to_netcdf(save_file, encoding=encoding, format='netCDF4', engine='netcdf4', unlimited_dims='time')


if __name__ == '__main__':
    # 'rucool/Saba/gliderdata/ru_ph_ncei_datasets' 'rucool/Saba/gliderdata/sbu/from_DAC' 'rucool/Saba/gliderdata/umaine/from_DAC'
    files = 'rucool/Saba/gliderdata/ru_ph_ncei_datasets/2022'
    save_directory = 'rucool/Saba/NOAA_SOE/data/output_nc/gliders'
    main(files, save_directory)
