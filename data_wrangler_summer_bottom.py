#!/usr/bin/env python

"""
Author: Lori Garzio on 10/18/2021
Last modified: 11/16/2022
Grab summer bottom-water pH and omega data from CODAP-NA, EcoMon, and glider datasets.
CODAP-NA dataset documented here: https://essd.copernicus.org/articles/13/2777/2021/
"""

import os
import numpy as np
import pandas as pd
import xarray as xr
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point
import statistics
import PyCO2SYS as pyco2
import pickle
import functions.common as cf
pd.set_option('display.width', 320, "display.max_columns", 20)  # for display in pycharm console


def main(lon_bounds, lat_bounds, codap_file, ecomon_files, glider_files, savedir):
    bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GEBCO_2014_2D_-100.0_0.0_-10.0_50.0.nc'
    extent = [-78, -65, 35, 45]
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0] - .1, extent[1] + .1),
                      lat=slice(extent[2] - .1, extent[3] + .1))

    # get CODAP data
    ds = xr.open_dataset(codap_file)
    idx = []
    codap_vars = dict(Month_UTC=np.array([]),
                      Year_UTC=np.array([]),
                      Cruise_ID=np.array([]),
                      Profile_number=np.array([]),
                      Latitude=np.array([]),
                      Longitude=np.array([]),
                      CTDPRES=np.array([]),
                      Depth=np.array([]),
                      Depth_bottom=np.array([]),
                      CTDTEMP_ITS90=np.array([]),
                      CTDTEMP_flag=np.array([]),
                      recommended_Salinity_PSS78=np.array([]),
                      recommended_Salinity_flag=np.array([]),
                      pH_TS_insitu_calculated=np.array([]),
                      pH_TS_insitu_measured=np.array([]),
                      pH_flag=np.array([]),
                      TALK=np.array([]),
                      TALK_flag=np.array([]),
                      Aragonite=np.array([]))

    # make sure the data are within the defined extent
    for i, lon in enumerate(ds.Longitude.values):
        if Polygon(list(zip(lon_bounds, lat_bounds))).contains(Point(lon, ds.Latitude.values[i])):
            idx.append(i)
            for key in codap_vars.keys():
                if key == 'Cruise_ID':
                    cid = ds.Cruise_ID.values[:, i]
                    cid = [x.decode('UTF-8') for x in cid]
                    codap_vars[key] = np.append(codap_vars[key], ''.join(cid).strip())
                else:
                    codap_vars[key] = np.append(codap_vars[key], ds[key].values[i])

    # select data from June - Aug
    df = pd.DataFrame(codap_vars)
    df = df[df.Month_UTC > 5]
    df = df[df.Month_UTC < 9]
    df['pH'] = df['pH_TS_insitu_measured']

    # use calculated pH if measured isn't available
    for idx, row in df.iterrows():
        if row['pH'] == -999:
            df.loc[row.name, 'pH'] = row.pH_TS_insitu_calculated

    df = df[df.pH != -999]

    # get rid of questionable pH flags
    df = df[df.pH_flag != 3]

    df['aragonite_final'] = ''

    # If aragonite saturation state is available, use that value. If it's not available, calculate aragonite saturation
    # state
    for idx, row in df.iterrows():
        if row.Aragonite == -999:
            # WOCE flags: 2 = Acceptable, 3 = Questionable, 6 = Average of duplicates, 9 = Missing
            if np.logical_or(row.TALK_flag == 2, row.TALK_flag == 6):  # if TA flags are acceptable
                if row.pH_flag != 3:  # if pH flags are acceptable
                    omega_arag, pco2, revelle = cf.run_co2sys_ta_ph(row.TALK,
                                                                    row.pH,
                                                                    row.recommended_Salinity_PSS78,
                                                                    row.CTDTEMP_ITS90,
                                                                    row.CTDPRES)
                    df.loc[row.name, 'aragonite_final'] = omega_arag
                else:
                    df.loc[row.name, 'aragonite_final'] = np.nan
            else:
                df.loc[row.name, 'aragonite_final'] = np.nan
        else:
            df.loc[row.name, 'aragonite_final'] = row.Aragonite

    # initialize dictionary to append bottom pH and aragonite data from cruises
    data = dict(cruise=np.array([]),
                depth=np.array([]),
                lat=np.array([]),
                lon=np.array([]),
                pH=np.array([]),
                aragonite=np.array([]))

    # for each profile on each cruise, find the deepest sample and compare the sample depth to either the water depth
    # provided in the file, or depth from the global bathymetry file
    cruises = np.unique(df.Cruise_ID)
    final_cruises = []
    for cruise in cruises:
        dfc = df[df.Cruise_ID == cruise]
        profile_num = np.unique(dfc.Profile_number)
        for profile in profile_num:
            dfc_profile = dfc[dfc.Profile_number == profile]
            maxdepth = np.nanmax(dfc_profile.Depth)
            if maxdepth > 10:  # sampling depth has to be >10m
                dfc_profile_max = dfc_profile[dfc_profile.Depth == maxdepth]
                profile_coords = [dfc_profile_max.Longitude.values[0], dfc_profile_max.Latitude.values[0]]
                if dfc_profile_max.Depth_bottom.values[0] != -999:  # compare to the recorded station depth
                    station_water_depth = dfc_profile_max.Depth_bottom.values[0]
                else:  # compare to the global bathymetry file
                    lat_idx = abs(bathy.lat.values - profile_coords[1]).argmin()
                    lon_idx = abs(bathy.lon.values - profile_coords[0]).argmin()
                    station_water_depth = -bathy.elevation[lat_idx, lon_idx].values

                # if the sample is within +/- 20% of the water column, keep the value
                depth_threshold = [station_water_depth * .8, station_water_depth * 1.2]
                if np.logical_and(maxdepth > depth_threshold[0], maxdepth < depth_threshold[1]):
                    data['cruise'] = np.append(data['cruise'], cruise)
                    data['depth'] = np.append(data['depth'], maxdepth)
                    data['lat'] = np.append(data['lat'], profile_coords[1])
                    data['lon'] = np.append(data['lon'], profile_coords[0])
                    data['pH'] = np.append(data['pH'], dfc_profile_max.pH.values[0])
                    data['aragonite'] = np.append(data['aragonite'], dfc_profile_max.aragonite_final.values[0])
                    final_cruises.append(cruise)

    # additional EcoMon data that aren't included in CODAP
    for ef in ecomon_files:
        df_ecomon = pd.read_csv(ef)
        df_ecomon = df_ecomon[df_ecomon.Month_UTC > 5]
        df_ecomon = df_ecomon[df_ecomon.Month_UTC < 9]
        df_ecomon.replace(-999, np.nan, inplace=True)
        df_ecomon.dropna(subset=['pH_TS_20C'], inplace=True)

        # the only time bottom depth isn't recorded is in the flow-thru data, so remove those rows
        df_ecomon.dropna(subset=['Depth_Bottom_meters'], inplace=True)

        # combine Station_ID and Cast_number to get a unique profile
        df_ecomon['profile_id'] = df_ecomon.Station_ID.astype(str) + '_' + df_ecomon.Cast_number.astype(str)

        # for each profile, find the deepest sample and compare the sample depth to the water column depth
        # keep the sample if it's within 20% of the bottom depth
        profile_num = np.unique(df_ecomon.profile_id)
        for profile in profile_num:
            df_ecomon_profile = df_ecomon[df_ecomon.profile_id == profile]
            maxdepth = np.nanmax(df_ecomon_profile.Depth_meters)
            if maxdepth > 10:  # sampling depth has to be >10m
                df_ecomon_profile_max = df_ecomon_profile[df_ecomon_profile.Depth_meters == maxdepth]
                profile_coords = [df_ecomon_profile_max.Longitude_Dec_Deg.values[0],
                                  df_ecomon_profile_max.Latitude_Dec_Deg.values[0]]
                station_water_depth = np.nanmax(df_ecomon_profile_max.Depth_Bottom_meters)
                depth_threshold = [station_water_depth * .8, station_water_depth * 1.2]
                if np.logical_and(maxdepth > depth_threshold[0], maxdepth < depth_threshold[1]):
                    data['cruise'] = np.append(data['cruise'], df_ecomon_profile_max.Cruise_ID.values[0])
                    data['depth'] = np.append(data['depth'], maxdepth)
                    data['lat'] = np.append(data['lat'], profile_coords[1])
                    data['lon'] = np.append(data['lon'], profile_coords[0])

                    # calculate pH at in situ temperature, pressure, and salinity using PyCO2SYS
                    ph_20c = df_ecomon_profile_max.pH_TS_20C.values[0]
                    par1_type = 3
                    ta = df_ecomon_profile_max['TA_umol/kg'].values[0]
                    if ta == -999:
                        ta = 2200
                    par2_type = 1

                    kwargs = dict(salinity=df_ecomon_profile_max.CTDSAL_PSS78.values[0],
                                  temperature=20,
                                  temperature_out=df_ecomon_profile_max.CTDTEMP_ITS90.values[0],
                                  pressure=0,
                                  pressure_out=df_ecomon_profile_max.CTDPRES_dbar.values[0],
                                  opt_pH_scale=1,
                                  opt_k_carbonic=4,
                                  opt_k_bisulfate=1,
                                  opt_total_borate=1,
                                  opt_k_fluoride=2)

                    results = pyco2.sys(ph_20c, ta, par1_type, par2_type, **kwargs)
                    data['pH'] = np.append(data['pH'], results['pH_out'])

                    # calculate aragonite
                    omega_arag, pco2, revelle = cf.run_co2sys_ta_ph(df_ecomon_profile_max['TA_umol/kg'].values[0],
                                                                    results['pH_out'],
                                                                    df_ecomon_profile_max.CTDSAL_PSS78.values[0],
                                                                    df_ecomon_profile_max.CTDTEMP_ITS90.values[0],
                                                                    df_ecomon_profile_max.CTDPRES_dbar.values[0])

                    data['aragonite'] = np.append(data['aragonite'], omega_arag)

    pickle_name = os.path.join(savedir, 'cruise_summer_bottom_data.pickle')
    with open(pickle_name, 'wb') as handle:
        pickle.dump(data, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # get glider data
    # for each glider, grab the data in the bottom 1m of each profile except when the profile is deeper than a
    # defined threshold (when the glider is off the shelf and no longer sampling the bottom)

    max_depth_threshold = dict(sbu01=170,
                               ru30=170,
                               um_242=300)

    glider_data = dict(pressure=np.array([]),
                       depth=np.array([]),
                       lat=np.array([]),
                       lon=np.array([]),
                       pH=np.array([]),
                       aragonite=np.array([]))

    for gf in glider_files:
        glider = gf.split('/')[-1].split('-')[0]
        ds = xr.open_dataset(gf)
        try:
            profileid = np.unique(ds.profile_time.values)
            profile_identifier = 'profile_time'
            phvar = 'ph_total_shifted'
        except AttributeError:
            profileid = np.unique(ds.profile_id.values)
            profile_identifier = 'profile_id'
            phvar = 'pH'
        try:
            arag = ds.saturation_aragonite.values
        except AttributeError:
            arag = ds.aragonite_saturation_state.values
        for pid in profileid:
            pidx = np.where(ds[profile_identifier].values == pid)[0]
            maxpress = np.nanmax(ds.pressure.values[pidx])
            if np.logical_and(maxpress > 10, maxpress < max_depth_threshold[glider]):
                idx = np.where(ds.pressure.values[pidx] > maxpress - 1)[0]

                # compare the glider depth to the global bathymetry file
                gl_lat = np.nanmean(ds.latitude.values[pidx][idx])
                gl_lon = np.nanmean(ds.longitude.values[pidx][idx])
                profile_coords = [gl_lon, gl_lat]
                lat_idx = abs(bathy.lat.values - profile_coords[1]).argmin()
                lon_idx = abs(bathy.lon.values - profile_coords[0]).argmin()
                station_water_depth = -bathy.elevation[lat_idx, lon_idx].values

                # if the glider sample is within +/- 20% of the water column, keep the value
                depth_threshold = [station_water_depth * .8, station_water_depth * 1.2]
                if np.logical_and(maxpress > depth_threshold[0], maxpress < depth_threshold[1]):
                    if np.sum(~np.isnan(arag[pidx][idx])):
                        glider_data['pressure'] = np.append(glider_data['pressure'], maxpress)
                        glider_data['depth'] = np.append(glider_data['depth'], np.nanmax(ds.depth.values[pidx][idx]))
                        glider_data['lat'] = np.append(glider_data['lat'], gl_lat)
                        glider_data['lon'] = np.append(glider_data['lon'], gl_lon)
                        glider_data['pH'] = np.append(glider_data['pH'], statistics.median(ds[phvar].values[pidx][idx]))
                        glider_data['aragonite'] = np.append(glider_data['aragonite'], statistics.median(arag[pidx][idx]))

    pickle_name = os.path.join(savedir, 'glider_summer_bottom_data.pickle')
    with open(pickle_name, 'wb') as handle:
        pickle.dump(glider_data, handle, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    lons = [-78, -65, -65, -78]
    lats = [35, 35, 45, 45]
    codap = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/CODAP_NA_v2021.nc'
    ecomon = ['/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/EcoMon/33GG20190815-GU1902_data.csv',
              '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/EcoMon/33HH20190522-HB1902_data.csv',
              '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/EcoMon/334B20210805-PC2104_Data.csv']
    gliders = ['/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/nc_files/sbu01-20210720T1628-profile-sci-delayed-qc_shifted_co2sys_final.nc',
               '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/nc_files/ru30-20190717T1812-profile-sci-delayed-dac.nc',
               '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/nc_files/ru30-20210716T1804-profile-sci-delayed_qc_final.nc',
               '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/nc_files/um_242-20210630T1916-profile-sci-delayed_shifted_final.nc']
    save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/pickle'
    main(lons, lats, codap, ecomon, gliders, save_directory)
