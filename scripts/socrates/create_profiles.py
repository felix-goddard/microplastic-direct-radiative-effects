import re
import sys
import scipy
import netCDF4
import numpy as np
import pandas as pd
import xarray as xr
from pathlib import Path

# Add the lib folder to the path so we can import constants
base_path = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(base_path))

import lib.constants as constants

profiles_dir = base_path / 'data' / 'socrates' / 'atmospheric_profiles'


def write_dataset(data, units, title, filename):        
    # Write a dataset to a netCDF file for use in SOCRATES

    print(f'Writing {filename}.nc')

    ncdf_file = netCDF4.Dataset(profiles_dir / f'{filename}.nc',
                                'w', format='NETCDF3_CLASSIC', clobber='true')

    def create_dim(vals, name, vtype, dims, units, title):
        dimension = ncdf_file.createDimension(name, len(vals))
        variable = ncdf_file.createVariable(name, vtype, dims)
        variable.units = units
        variable.title = title
        variable[:] = vals

    def create_var(vals, name, vtype, dims, units, title):
        variable = ncdf_file.createVariable(name, vtype, dims)
        variable.units = units
        variable.title = title
        variable[:] = vals

    create_dim(data.lon.values, 'lon', 'f4', 'lon', 'degree', 'LONGITUDE')
    create_dim(data.lat.values, 'lat', 'f4', 'lat', 'degree', 'LATITUDE')

    if 'plev' in data.dims:
        data = data.transpose('plev', 'lon', 'lat')
        data_dims = ['plev', 'lon', 'lat']
        create_dim(data.plev.values, 'plev', 'f4', 'plev', 'Pa', 'PRESSURE')
    elif 'basis' in data.dims:
        data = data.transpose('basis', 'lon', 'lat')
        data_dims = ['basis', 'lon', 'lat']
        create_dim(data.basis.values, 'basis', 'i2', 'basis', 'none', 'BASIS FUNCTION')
    else:
        data = data.transpose('lon', 'lat')
        data_dims = ['lon', 'lat']
    
    create_var(data.values, data.name, 'f4', data_dims, units, title)

    ncdf_file.close()


#=================================================================================================================
# Create profiles of atmospheric and radiative conditions
#=================================================================================================================

atmos_src = xr.open_dataset(base_path / 'data' / 'socrates' / 'atmos_climatology.nc') # 3D atmospheric data
surf_src = xr.open_dataset(base_path / 'data' / 'socrates' / 'surface_climatology.nc') # 2D surface data
solar_irradiance = xr.open_dataset(base_path / 'data' / 'socrates' / 'solar_irradiance.nc').stoa


def interp_linear_logp(var):
    # Interpolate a variable linearly in log(pressure) space

    mid_pressures = var.plev.rolling(plev=2).mean().dropna('plev')
    interp = scipy.interpolate.interp1d(
        np.log(var.plev.values), var.values,
        axis=var.dims.index('plev'))(np.log(mid_pressures.values))
    return xr.DataArray(
        data=interp, dims=var.dims,
        coords=dict(
            plev=('plev', mid_pressures.values),
            lon=('lon', var.lon.values),
            lat=('lat', var.lat.values),
        ),
    ).rename(var.name)


def interp_log_logp(var):
    # Interpolate the log of a variable linearly in log(pressure) space

    mid_pressures = var.plev.rolling(plev=2).mean().dropna('plev')
    interp = np.exp(scipy.interpolate.interp1d(
        np.log(var.plev.values), np.log(var.values),
        axis=var.dims.index('plev'))(np.log(mid_pressures.values)))
    return xr.DataArray(
        data=interp, dims=var.dims,
        coords=dict(
            plev=('plev', mid_pressures.values),
            lon=('lon', var.lon.values),
            lat=('lat', var.lat.values),
        ),
    ).rename(var.name)


# Calculate atmospheric fields ----------------------------------------------------------------------------------

t_centers = interp_linear_logp(atmos_src.t)
write_dataset(atmos_src.t, 'K', 'TEMPERATURE', 'temperature_layers')
write_dataset(t_centers,   'K', 'TEMPERATURE', 'temperature_centers')

write_dataset(interp_log_logp(atmos_src.q),  'None', 'MMR OF WATER VAPOUR', 'water_vapour')
write_dataset(interp_log_logp(atmos_src.o3), 'None', 'MMR OF OZONE', 'o3')

atmospheric_field = lambda name, value: xr.full_like(t_centers, value).rename(name)
write_dataset(atmospheric_field('o2', 0.2314),    'None', 'O2 MMR', 'o2')
write_dataset(atmospheric_field('co2', 5.241e-4), 'None', 'CO2 MMR', 'co2')

# Calculate surface fields --------------------------------------------------------------------------------------

write_dataset(atmos_src.t.sel(plev=100000).rename('tstar').expand_dims(plev=[100000]),
              'K', 'SURFACE TEMPERATURE', 'surface_temperature')

albedo = lambda x: x.expand_dims(basis=[1]).rename('alb')
write_dataset(albedo(surf_src.fal),                  'None', 'SURFACE ALBEDO', 'albedo_sw')
write_dataset(xr.full_like(albedo(surf_src.fal), 0), 'None', 'SURFACE ALBEDO', 'albedo_lw')


# Calculate solar zenith angles and solar irradiance ------------------------------------------------------------

# Only do this if the files don't already exist, since this is slow
if not ((profiles_dir / 'solar_zenith.nc').exists() and (profiles_dir / 'solar_irradiance.nc').exists()):

    time = np.arange(0, 365, .5/24)
    lon, lat, time = np.meshgrid(surf_src.lon.values, surf_src.lat.values, time)
    
    day_number = np.floor(time)
    N = 2*np.pi * day_number / 365
    
    equation_of_time = (7.5e-5 + 1.868e-3 * np.cos(N) - 3.2077e-2 * np.sin(N)
                        - 1.4615e-2 * np.cos(2*N) - 4.0849e-2 * np.sin(2*N))
    
    hour_angle = ((time-day_number)*(24/12) - 1 + lon/180) * np.pi + equation_of_time
    
    solar_declination = (6.918e-3 - 0.399912*np.cos(N) + 0.070257*np.sin(N)
                         - 6.758e-3*np.cos(2*N) + 9.07e-4*np.sin(2*N)
                         - 2.697e-3*np.cos(3*N) + 1.480e-3*np.sin(3*N))
    
    cos_solar_zenith_angle = (np.sin(solar_declination) * np.sin(np.deg2rad(lat))
                              + np.cos(solar_declination) * np.cos(np.deg2rad(lat)) * np.cos(hour_angle))
    
    avg = np.average(cos_solar_zenith_angle, axis=2, weights=cos_solar_zenith_angle>0)
    avg = np.rad2deg(np.arccos(np.average(avg, axis=1)))
    
    szen = xr.full_like(solar_irradiance, 0).rename('szen')
    szen[:] = np.tile(avg, (len(surf_src.lon), 1)).T

    write_dataset(szen, 'degree', 'SOLAR ZENITH', 'solar_zenith')

    average_irradiance = 340 # Wm-2
    scaling_factor = average_irradiance / ( solar_irradiance*np.cos(np.deg2rad(szen)) ).weighted(np.cos(np.deg2rad(szen.lat))).mean()
    write_dataset((scaling_factor*solar_irradiance).rename('stoa'), 'Wm-2', 'SOLAR IRRADIANCE', 'solar_irradiance')


#=================================================================================================================
# Create profiles of MnP particle concentrations
#=================================================================================================================

layer_t = (xr.open_dataset(profiles_dir / 'temperature_layers.nc').t.transpose('plev','lat','lon'))
center_t = (xr.open_dataset(profiles_dir / 'temperature_centers.nc').t.transpose('plev','lat','lon'))
center_q = (xr.open_dataset(profiles_dir / 'water_vapour.nc').q.transpose('plev','lat','lon'))

virtual_temperature = center_t * (1 + center_q * (1/constants.water_vapour_mass_ratio - 1))
log_pressure_diff = np.log(layer_t.plev).diff('plev').assign_coords(plev=layer_t.plev-layer_t.plev.diff('plev')/2)
thickness = (constants.dry_gas_constant/constants.gravity) * virtual_temperature * log_pressure_diff
zl = np.cumsum(np.insert(thickness.values[::-1,:,:], [0], 0, axis=0), axis=0)[::-1]
z = xr.zeros_like(center_t)
z.values[:] = (zl[1:,:,:] + zl[:-1,:,:])/2


def vertically_exponential(top_height, base=0.3):
    fraction = base**(z/10e3)
    fraction.values[z.values > top_height] = 0
    return fraction


def vertically_wellmixedBL():
    bl_height = 2e3
    tropopause_height = (
        xr.open_dataset(base_path / 'data' / 'tropopause_height.nc')
        .tropopause_height.mean('time').transpose('lat','lon').values)
    fraction = np.minimum(1, 0.3**((z - bl_height)/10e3))
    fraction.values[z.values > tropopause_height] = 0
    return fraction


def horizontally_uniform(conc):
    return conc # MPs/m3


def horizontally_evangeliou(conc):
    mp_conc = xr.open_dataset(base_path / 'data' / 'nonuniform_mp_concentration.nc').mp_concentration
    mp_conc = conc * mp_conc / mp_conc.weighted(np.cos(np.deg2rad(mp_conc.lat))).mean()
    return mp_conc



template = xr.zeros_like(center_t).rename('mp_conc')

for conc in [1]:
    template.values[:] = vertically_exponential(2e3) * horizontally_uniform(conc)
    write_dataset(template, 'M-3', 'MP CONCENTRATION', f'mps_exponential(0.3,2km)_uniform({conc}m-3)')
    
    template.values[:] = vertically_wellmixedBL() * horizontally_uniform(conc)
    write_dataset(template, 'M-3', 'MP CONCENTRATION', f'mps_wellmixedBL_uniform({conc}m-3)')

    template.values[:] = vertically_exponential(2e3) * horizontally_evangeliou(conc)
    write_dataset(template, 'M-3', 'MP CONCENTRATION', f'mps_exponential(0.3,2km)_evangeliou({conc}m-3)')
    
    template.values[:] = vertically_wellmixedBL() * horizontally_evangeliou(conc)
    write_dataset(template, 'M-3', 'MP CONCENTRATION', f'mps_wellmixedBL_evangeliou({conc}m-3)')
