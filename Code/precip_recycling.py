#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
;;#############################################################################
;;
;; precipitation_recycling.py
;; Author: Matteo Mastropierro (mastropierromatteo@gmail.com)
;; Ca' Foscari University of Venice, ITALY
;;
;;#############################################################################
;;
;; Description
;;    This script is used to calculate the Budyko recycling coefficient (B).
;;    For a given region, B is the ratio of locally evaporated to advected
;;    moisture.
;;
;;    Using the equation from Brukbaker et al. (1993):
;;
;;    B = 1 + EA/2Fin
;;
;;    where E is evap per unit area, A is the area of the region, and Fin is
;;    the total moisture flux into the region.
;;
;;    B can be used to estimate the precipitation recycling ratio (Rr), i.e. the
;;    proportion of precipitation within an area that originated from
;;    evaporation within that area.
;;
;;    Rr = 1 - (1/B)
;;
;; Requirements
;;
;;#############################################################################
"""

import os
import glob
import math
import numpy as np
import xarray as xr
import pickle

def find_dl(lat, lons):
    """
    Find average length of 1 grid cell - convert lat to radians, take cos +
    divide by Earth's circumference (m)
    """
    earth_circ = (40075*10**3)  # Earth's circumference in m
    lon_shape = len([lons])
    lat_rad = np.radians(abs(lat))
    dl = math.cos(lat_rad) * (earth_circ / lon_shape)  # length of pixel in m
    return dl

def lon180(ds):
    ds.coords['lon'] = (ds.coords['lon'] + 180) % 360 - 180
    ds = ds.sortby(ds.lon)
    return ds

def get_pixel_size(lat, lon):
    if lat[0] > lat[-1]:
        temp_lat = lat[::-1]
    else:
        temp_lat = lat
    r = 6.371 * 1e6
    rad = (2 * math.pi / 360)  # (m)
    da = np.nan * np.zeros((len(temp_lat)))  # (m^2)

    for i in range(len(temp_lat) - 1):
        da[i] = (2 * math.pi * (1 / len(lon)) *
                 r ** 2 * (math.sin(rad * temp_lat[i + 1]) -
                           math.sin(rad * temp_lat[i])))

    # Check if top and bottom latitude are same
    if temp_lat[0] == -temp_lat[-1]:
        da[-1] = da[0]
    return da


def main():
    data_path = "F:/Data/LUMIP"
    models = ["UKESM1-0-LL"]

    for i, mm in enumerate(models):
        
        print("Importing data of " + mm)
        # Open data
        q1 = xr.open_dataset(glob.glob(os.path.join(data_path, mm, f"hus_Amon_{mm}_ssp370_*.nc"))[0], engine='netcdf4', chunks={"time": 258})["hus"]
        u1 = xr.open_dataset(glob.glob(os.path.join(data_path, mm, f"ua_Amon_{mm}_ssp370_*.nc"))[0], engine='netcdf4', chunks={"time": 258})["ua"]
        v1 = xr.open_dataset(glob.glob(os.path.join(data_path, mm, f"va_Amon_{mm}_ssp370_*.nc"))[0], engine='netcdf4', chunks={"time": 258})["va"]

        q2 = xr.open_dataset(glob.glob(os.path.join(data_path, mm, f"hus_Amon_{mm}_ssp370-ssp126Lu_*.nc"))[0], engine='netcdf4', chunks={"time": 258})["hus"]
        u2 = xr.open_dataset(glob.glob(os.path.join(data_path, mm, f"ua_Amon_{mm}_ssp370-ssp126Lu_*.nc"))[0], engine='netcdf4', chunks={"time": 258})["ua"]
        v2 = xr.open_dataset(glob.glob(os.path.join(data_path, mm, f"va_Amon_{mm}_ssp370-ssp126Lu_*.nc"))[0], engine='netcdf4', chunks={"time": 258})["va"]

        data_path = "G:/My Drive/MPIM/data/"

        scenario = 'ssp126Lu'
        filepath = glob.glob(os.path.join(data_path, f'xr_{scenario}_{mm}_pft.nc'))[0]
        xr_aff_pft = xr.open_dataset(filepath, drop_variables=["time_bnds", "lon_bnds", "lat_bnds"], engine='netcdf4', chunks={"time": 258})["treeFrac"]
        filepath = glob.glob(os.path.join(data_path, f'xr_{scenario}_{mm}.nc'))[0]
        xr_aff = xr.open_dataset(filepath, drop_variables=["time_bnds", "lon_bnds", "lat_bnds"], engine='netcdf4', chunks={"time": 258})["evspsbl"]

        scenario = 'ssp370'
        filepath = glob.glob(os.path.join(data_path, f'xr_{scenario}_{mm}_pft.nc'))[0]
        xr_ctl_pft = xr.open_dataset(filepath, drop_variables=["time_bnds", "lon_bnds", "lat_bnds"], engine='netcdf4', chunks={"time": 258})["treeFrac"]
        filepath = glob.glob(os.path.join(data_path, f'xr_{scenario}_{mm}.nc'))[0]
        xr_ctl = xr.open_dataset(filepath, drop_variables=["time_bnds", "lon_bnds", "lat_bnds"], engine='netcdf4', chunks={"time": 258})["evspsbl"]

        # Preprocess data
        xr_delta_pft = (xr_aff_pft.sel(time=slice("2071-01", "2100-12")).mean(dim="time") - xr_ctl_pft.sel(time=slice("2071-01", "2100-12")).mean(dim="time"))
        treefrac_pos = xr_delta_pft.where(xr_delta_pft > 10)

        et = (xr_aff - xr_ctl)  # Evapotranspiration in Kg m-2 s-1
        q2.assign_coords(time = q2.time); q2 = q2.assign_coords(plev=q1.plev)
        u2.assign_coords(time = u2.time); u2 = u2.assign_coords(plev=u1.plev)
        v2.assign_coords(time = v2.time); v2 = v2.assign_coords(plev=v1.plev)

        q = (q2 - q1)
        u = (u2 - u1)
        v = (v2 - v1)

        q = q.sel(lat=slice(-61,90))
        u = u.sel(lat=slice(-61,90))
        v = v.sel(lat=slice(-61,90))
        if treefrac_pos.lat.shape != q.lat.shape:       
            q = q.sel(lat=slice(-60,90))
            u = u.sel(lat=slice(-60,90))
            v = v.sel(lat=slice(-60,90))

        treefrac_pos = treefrac_pos.assign_coords(lat=q.lat); treefrac_pos = treefrac_pos.assign_coords(lon=q.lon)
        et = et.assign_coords(lat=q.lat); et = et.assign_coords(lon=q.lon); et = et.assign_coords(time=q.time)

        print("Computing vertical flux integration")
        # Flux vertical integration
        # Ensure that the datasets have matching dimensions
        assert q.dims == u.dims == v.dims

        # Pressure levels (Pa)
        pressure = q.coords['plev']

        # Calculate the zonal (qx) and meridional (qy) water vapor fluxes (kg/(m^2*s))
        qx = q * u
        qy = q * v

        # Integrate vertically over all pressure levels
        dp = np.diff(pressure)  # Difference between consecutive pressure levels

        # Add an extra level at the top for integration purposes
        dp = np.append(dp, dp[-1])

        # Reshape dp to make it broadcastable across qx and qy
        dp = dp.reshape((1, -1, 1, 1))

        # Vertical integration
        qu = (qx * dp).sum(dim='plev') / 9.81  # Dividing by gravitational acceleration to convert to (kg/m^2/s)
        qv = (qy * dp).sum(dim='plev') / 9.81

        # Resulting vertically integrated water vapor fluxes
        print("Zonally integrated fluxes qu:  ", qu)
        print("Meridionally integrated fluxes qv: ", qv)

        pixel_size_grid = get_pixel_size(qu.lat.values, qu.lon.values)
        pixel_size_grid = np.array([pixel_size_grid] * len(qu.lon)).transpose()

        qu = qu.compute()
        qv = qv.compute()
        et = et.compute()

        # Compute the operation on the cells with at least 10% of increased treeFrac in the 2070-2100 period
        cells = treefrac_pos.stack(cell=["lon", "lat"]).where(treefrac_pos.stack(cell=["lon", "lat"]).notnull(), drop=True).cell
        grids = cells.values
        res_y = (treefrac_pos.lat[1] - treefrac_pos.lat[0]).values
        res_x = (treefrac_pos.lon[1] - treefrac_pos.lon[0]).values

        qu_stack = qu.stack(cell=["lon", "lat"])
        qv_stack = qv.stack(cell=["lon", "lat"])
        et_stack = et.stack(cell=["lon", "lat"])

        xr_B = xr.DataArray(data=None, coords=[qu_stack.time, qu_stack.cell], dims=["time", "cell"])
        xr_B["B"] = xr.DataArray(data=None, coords=[qu_stack.time, qu_stack.cell], dims=["time", "cell"])
        xr_B["Rr"] = xr.DataArray(data=None, coords=[qu_stack.time, qu_stack.cell], dims=["time", "cell"])

        print("Computing vPrecipitation recycling raation as from Brubaker")
        # Iterate Precipitation recycling ratio over cells of Delta_treefrac>10%
        for i, c in enumerate(cells):
            locator = {'cell': c}
            g = grids[i]

            # Lon & Lat indices
            ilon = np.where(qu.lon == g[0])[0][0]
            ilat = np.where(qu.lat == g[1])[0][0]

            lonmin = g[0] - res_x / 2
            lonmax = g[0] + res_x / 2
            latmin = g[1] - res_y / 2
            latmax = g[1] + res_y / 2

            # Find average grid cell length on E, W, N and S sides of domain in m
            dl_e = find_dl(np.mean((latmin, latmax)), lonmax)
            dl_w = find_dl(np.mean((latmin, latmax)), lonmin)
            dl_n = find_dl(latmax, np.mean((lonmin, lonmax)))
            dl_s = find_dl(latmin, np.mean((lonmin, lonmax)))

            # Calculate inward moisture flux in kg s-1
            inflow_from_E = -dl_e * qu_stack.loc[locator]
            inflow_from_W = dl_w * qu_stack.loc[locator]
            inflow_from_N = -dl_n * qv_stack.loc[locator]
            inflow_from_S = dl_s * qv_stack.loc[locator]

            # Initialize flux_in and flux_out with zeros
            flux_in = xr.zeros_like(inflow_from_E)
            flux_out = xr.zeros_like(inflow_from_E)

            # Calculate total water vapour flux into the cell
            flux_in += inflow_from_E.where(inflow_from_E > 0, 0)
            flux_in += inflow_from_W.where(inflow_from_W > 0, 0)
            flux_in += inflow_from_N.where(inflow_from_N > 0, 0)
            flux_in += inflow_from_S.where(inflow_from_S > 0, 0)

            # Calculate Budyko recycling coefficient
            A = pixel_size_grid[ilat, ilon]
            EA = et_stack.loc[locator] * A

            B = 1 + (EA.values / (2 * flux_in.values))
            Rr = 1 - (1 / B)

            xr_B["B"].loc[locator] = xr.DataArray(data=B)
            xr_B["Rr"].loc[locator] = xr.DataArray(data=Rr)

        xr_B = xr_B.unstack()
        xr_B["B"] = xr_B["B"].astype(np.float64)
        xr_B["Rr"] = xr_B["Rr"].astype(np.float64)

        # Save and export regression list data
        print("Saving to: ", data_path)
        with open(os.path.join(data_path, f"budyko_pr_{mm}"), "wb") as fp:
            pickle.dump(xr_B, fp)

        # Free up memory
        del(q,u,v,q1,q2,u1,u2,v1,v2,et)

if __name__ == "__main__":
    main()

