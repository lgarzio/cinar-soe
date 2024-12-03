#!/usr/bin/env python

"""
Author: Lori Garzio on 11/29/2021
Last modified: 11/29/2021
Plot 3D maps of glider pH.
"""

import datetime as dt
import os
import xarray as xr
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # need this for 3D scatter plot
import matplotlib.pyplot as plt
from matplotlib import colors as c
import cmocean as cmo
plt.rcParams.update({'font.size': 14})


def main(save_dir, file):
    if len(file) == 1:  # UMaine deployment
        deployment = file[0].split('/')[-1].split('-profile-sci')[0]
        bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GMRTv3_9_20211116topo-gulfofmaine-lowres.grd'
        extent = [-71, -67, 41, 45]
        depth_max = 350
        yz_tickpad = 10
        yz_labelpad = 16
        rotate_view = [25, -80]
    else:  # ru30 and sbu01 on the same map
        deployment = 'ru30_sbu01'
        bathymetry = '/Users/garzio/Documents/rucool/bathymetry/GMRTv3_9_20211118topo-mab-highres.grd'
        extent = [-75, -71.5, 37.8, 41.2]
        depth_max = 250
        yz_tickpad = 6
        yz_labelpad = 12
        rotate_view = [25, -70]

    # manipulate bathymetry
    bathy = xr.open_dataset(bathymetry)
    bathy = bathy.sel(lon=slice(extent[0], extent[1]),
                      lat=slice(extent[2], extent[3]))
    x, y = np.meshgrid(bathy.lon.values, bathy.lat.values)
    alt = bathy.altitude
    lm_idx = np.logical_and(alt > 0, alt > 0)
    deep_idx = np.logical_and(alt < -depth_max, alt < -depth_max)

    # mask land elevation and deep values
    alt.values[lm_idx] = np.nan
    alt.values[deep_idx] = np.nan

    # find the coastline
    coast = np.empty(shape=np.shape(alt))
    coast[:] = np.nan

    # create a land variable
    landmask = alt.values.copy()
    landmask[lm_idx] = 0
    landmask[~lm_idx] = np.nan

    # create 3D plot
    fig, ax = plt.subplots(figsize=(12, 10), subplot_kw={"projection": "3d",
                                                         "computed_zorder": False})

    # change the scaling of the plot area
    ax.set_box_aspect(aspect=(2, 2, 1))

    # Remove gray panes and axis grid
    ax.xaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('white')
    ax.yaxis.pane.fill = False
    ax.yaxis.pane.set_edgecolor('white')
    ax.zaxis.pane.fill = False
    ax.zaxis.pane.set_edgecolor('white')
    ax.grid(False)

    # Create bins between cbar_min and cbar_max with a step of cbar_step
    cbar_min = 7.74
    cbar_max = 8.2
    cbar_step = .02
    ticks = np.arange(cbar_min, cbar_max + cbar_step, cbar_step)
    cmap = cmo.cm.matter

    # Use bins/ticks as boundarys for colormaps
    norm = c.BoundaryNorm(ticks, cmap.N)

    ax.plot_surface(x, y, -alt, color='lightgray', zorder=-1, alpha=.3)  # plot bathymetry
    ax.plot_wireframe(x, y, -alt, color='gray', zorder=0, alpha=.3)  # plot wireframe bathymetry
    ax.plot_surface(x, y, landmask, color='tan', zorder=1, shade=False)  # plot land

    for f in file:
        ds = xr.open_dataset(f)

        if 'sbu01' in f:
            start_date = dt.datetime(2021, 8, 14, 0, 0)
            end_date = dt.datetime(2021, 8, 21, 0, 0)
            ds = ds.sel(time=slice(start_date, end_date))

        try:
            sct = ax.scatter(ds.longitude.values, ds.latitude.values, ds.depth.values, c=ds.ph_total_shifted.values, s=2,
                             cmap=cmap, norm=norm, zorder=2)
        except AttributeError:
            sct = ax.scatter(ds.longitude.values, ds.latitude.values, ds.depth.values, c=ds.pH.values, s=2,
                             cmap=cmap, norm=norm, zorder=2)
    cbar = plt.colorbar(sct, label='pH', extend='both', shrink=0.7, pad=0.08)
    ax.invert_zaxis()
    ax.set_zlabel('Depth (m)', labelpad=yz_labelpad)
    ax.set_ylabel('Latitude', labelpad=yz_labelpad)
    ax.set_xlabel('Longitude', labelpad=14)

    # add space between ticks and tick labels
    ax.tick_params(axis='y', which='major', pad=yz_tickpad)
    ax.tick_params(axis='z', which='major', pad=yz_tickpad)

    if rotate_view:
        # ax.view_init(elev=30, azim=-60)  # defaults
        ax.view_init(elev=rotate_view[0], azim=rotate_view[1])  # rotate the view

    sfile = os.path.join(save_dir, f'{deployment}_summer2021_pH_3D-sbu01lastleg.png')
    plt.savefig(sfile, dpi=200)
    plt.close()


if __name__ == '__main__':
    save_directory = '/Users/garzio/Documents/rucool/Saba/NOAA_SOE2021/data/3d'
    #file = ['/Users/garzio/Documents/rucool/Saba/gliderdata/2021/um_242-20210630T1916/delayed/um_242-20210630T1916-profile-sci-delayed_shifted.nc']
    file = ['/Users/garzio/Documents/rucool/Saba/gliderdata/2021/ru30-20210716T1804/delayed/ru30-20210716T1804-profile-sci-delayed_shifted.nc',
            '/Users/garzio/Documents/rucool/Saba/gliderdata/2021/sbu01-20210720T1628/delayed/sbu01-20210720T1628-profile-sci-delayed_shifted.nc']
    main(save_directory, file)
