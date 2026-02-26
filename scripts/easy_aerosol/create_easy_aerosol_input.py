import re
import cftime
import numpy as np
import xarray as xr
from pathlib import Path

#==============================================================================================================

# Each of these entries defines a different combination of colour, size distribution,
#   and horizontal and vertical distributions. Directories will be created for each
#   to which the output files will be written.
experiments = [
    'clear_gamma(1-100)_exponential(0.3,2km)_uniform(420m-3)',  # Linearity test, target = 0.05 W/m2
    'clear_gamma(1-100)_exponential(0.3,2km)_uniform(840m-3)',  # Linearity test, target = 0.10 W/m2

    'clear_gamma(1-100)_exponential(0.3,2km)_uniform(1680m-3)',         # Experiment 1
    'clearUV_gamma(1-100)_exponential(0.3,2km)_uniform(1760m-3)',       # Experiment 2
    'black_gamma(1-100)_exponential(0.3,2km)_uniform(860m-3)',          # Experiment 3
    'mixedUV_gamma(1-100)_exponential(0.3,2km)_uniform(7230m-3)',       # Experiment 4
    'clear_powerlaw(1.52,1-100)_exponential(0.3,2km)_uniform(5710m-3)', # Experiment 5
    'clear_gamma(1-100)_wellmixedBL_uniform(2620m-3)',                  # Experiment 6
    'clear_gamma(1-100)_exponential(0.3,2km)_evangeliou(3070m-3)',      # Experiment 7
    'mixedUV_powerlaw(1.52,1-100)_wellmixedBL_evangeliou(650m-3)',      # Experiment 8
]

#==============================================================================================================

base_path = Path(__file__).resolve().parent.parent.parent

# This directory contains the files of optical properties
opt_prop_dir = base_path / 'data' / 'optical_properties'

# This is the netCDF file with the surface altitude data
orography_file = base_path / 'data' / 'easy_aerosol' / 'orography.nc'
orog = xr.open_dataset(orography_file).surface_altitude

# This file specifies the level structure of the model; this file
# does not exist in this repository, but can be created from the level
# set L85(50t,35s)85 given in namelist format at
#   https://gmd.copernicus.org/articles/10/1487/2017/gmd-10-1487-2017-supplement.pdf
levels_file = base_path / 'data' / 'easy_aerosol' / 'vertical_levels_file'

# These are the times we output at
time = cftime.date2num(
    [cftime.datetime(2000,m,16) for m in range(1,13)],
    'days since 2000-01-01 00:00:00', calendar='360_day'
)

#==============================================================================================================
# Handy helpful helper functions, for things like reading files and decoding coordinates


def retrieve_levels_data():
    with open(levels_file) as f:
        content = f.read()
        top_height = float(re.search(r'z_top_of_model\s*=\s*(\d+(:?\.\d+)?)\s*,', content).group(1))
        const_level_idx = int(re.search(r'first_constant_r_rho_level\s*=\s*(\d+(:?\.\d+)?)\s*,', content).group(1))

        eta_theta = re.search(r'eta_theta=\s*((?:(?:\d+(?:\.\d+)?(?:E[+-]\d+)?),\s*)+)', content, re.MULTILINE).group(1)
        eta_theta = [float(x.strip()) for x in eta_theta.split(',') if x.strip()]

        eta_rho = re.search(r'eta_rho=\s*((?:(?:\d+(?:\.\d+)?(?:E[+-]\d+)?),\s*)+)', content, re.MULTILINE).group(1)
        eta_rho = [float(x.strip()) for x in eta_rho.split(',') if x.strip()]

    # We subtract 1 from the constant rho level index because Fortran is 1-indexed
    return eta_rho, eta_rho[const_level_idx-1], top_height


def altitude_3d():
    eta, eta_ref, z_top = retrieve_levels_data()
    
    z_full = np.full((len(eta), *orog.shape), np.nan, np.float64)
    for i in range(len(eta)):
        # Calculation of altitude, based on the definition of the vertical coordinate given in the UMDP papers
        #   (https://code.metoffice.gov.uk/doc/um/latest/papers/umdp_F03.pdf)
        z_full[i,:,:] = eta[i]*z_top + orog.values * np.maximum(0, 1-eta[i]/eta_ref)**2

    return z_full # Return the height above sea level


def read_scattering_properties_file(filename): 
    with open(filename) as f:
        lines = f.readlines()
        while not lines[0].startswith('Band'):
            lines = lines[1:] # Skip all the preamble
        lines = lines[2:] # Skip the header lines
        lines = [x.strip().split()[1:] for x in lines if x.strip()]
        lines = [list(map(float, x)) for x in lines]
    
        col = lambda n: [x[n] for x in lines]
        absorption_coeff = col(0)
        scattering_coeff = col(1)
        asymmetry_factor = col(2)

    return absorption_coeff, scattering_coeff, asymmetry_factor


def make_3d_scattering_data(name, density, opt_prop_file):
    abs, sca, asy = read_scattering_properties_file(opt_prop_file)
    data_3d = np.zeros((len(time), len(asy), *density.shape[1:]))

    match name:
        case 'asymmetry':
            f = lambda i: asy[i]
        case 'extinction':
            f = lambda i: (abs[i] + sca[i]) * density
        case 'absorption':
            f = lambda i: abs[i] * density

    for i in range(len(asy)):
        data_3d[:,i,:] = f(i)

    return data_3d


def save_netcdf(path, data, var_name, units):
    ds = xr.Dataset(
        data_vars={
            var_name: (
                ['time', 'specband', 'model_level_number', 'latitude', 'longitude'],
                data,
                dict(vertical_scaling='all_levels',
                     standard_name=var_name,
                     long_name=var_name.replace('_', ' '),
                     units=units,))
        },
        coords=dict(
            time=(
                'time', time, dict(
                    calendar='360_day',
                    calendar_flexible=1,
                    standard_name='time',
                    units='days since 2000-01-01 00:00:00',)),
            specband=(
                'specband', range(1, 1+data.shape[1]), dict(
                    long_name='spectral band',
                    units='1',)),
            model_level_number=(
                ['model_level_number'], range(1, 1+data.shape[2]), dict(
                    long_name='model rho levels (Charney-Phillips grid)',
                    units='1',)),
            latitude=(
                ['latitude'], orog.latitude.values, dict(
                    standard_name='latitude',
                    units='degrees_north',)),
            longitude=(
                ['longitude'], orog.longitude.values, dict(
                    standard_name='longitude',
                    units='degrees_east',))
        ),
        attrs=dict(
            Conventions='CF-1.5',
            update_freq_in_hours='120',
            update_type='2',
        )
    )
    
    encoding = {var_name: {'_FillValue': np.nan}}
    ds.to_netcdf(path, encoding=encoding)


#==============================================================================================================
# Classes and functions that define the different assumptions we make


class VerticalDistribution():
    @staticmethod
    def exponential_scaling(top_height, base=0.3):
        altitude = altitude_3d()
        h = altitude - orog.values
        fraction = base**(h/10e3)
        fraction[altitude > top_height] = 0
        return fraction

    @staticmethod
    def well_mixed_boundary(bl_height):
        # The well-mixed boundary layer cuts off at the tropopause height, which is time-dependent, so we
        # must add an additional time axis to the distribution
        altitude = altitude_3d()
        tropopause_height = xr.open_dataset(base_path / 'data' / 'tropopause_height.nc').tropopause_height.values
        altitude = np.repeat(altitude[np.newaxis,:,:,:], tropopause_height.shape[0], axis=0)
        h = altitude - orog.values
        fraction = np.minimum(1, 0.3**((h - bl_height)/10e3))
        for i in range(altitude.shape[1]):
            fraction[:,i,:,:][altitude[:,i,:,:] > tropopause_height] = 0
        return fraction


class HorizontalDistribution():
    @staticmethod
    def constant(value):
        return np.ones_like(orog) * value

    @staticmethod
    def from_dataarray(arr, multiplier=1, lon_name='longitude', lat_name='latitude'):
        arr = arr.interp(**{lon_name: orog.longitude, lat_name: orog.latitude},
                         kwargs={'fill_value': np.nan})
        return arr.values


#==============================================================================================================


if __name__ == '__main__':
    for experiment in experiments:
        dir = base_path / 'data' / 'easy_aerosol' / 'forcing' / experiment
    
        colour, size, vertical, horizontal = experiment.split('_')
    
        # Optical properties --------------------------------------------------------------------------------------
    
        sw_file = Path(opt_prop_dir) / f'{colour}_{size}_sw'
        lw_file = Path(opt_prop_dir) / f'{colour}_{size}_lw'
    
        if not (sw_file.exists() and lw_file.exists()):
            print(f'No optical properties data exists for the combination of colour=`{colour}`',
                  f'and size=`{size}`; skipping `{experiment}`')
            continue
    
        # Vertical distribution -----------------------------------------------------------------------------------
    
        if (match := re.match(r'exponential\((\d+(?:\.\d+)?),(\d+(?:\.\d+)?)km\)', vertical)) is not None:
            vertical = VerticalDistribution.exponential_scaling(float(match.group(2)) * 1e3, base=float(match.group(1)))
    
        elif vertical == 'wellmixedBL':
            vertical = VerticalDistribution.well_mixed_boundary(2e3)
    
        else:
            print(f'No vertical distribution function exists for the option `{vertical}`; skipping `{experiment}`')
            continue
    
        # Horizontal distribution ---------------------------------------------------------------------------------
    
        if (match := re.match(r'uniform\((\d+(?:\.\d+)?)m-3\)', horizontal)) is not None:
            horizontal = HorizontalDistribution.constant(float(match.group(1)))
    
        elif (match := re.match(r'evangeliou\((\d+(?:\.\d+)?)m-3\)', horizontal)) is not None:
            target_density = float(match.group(1))
            mp_conc = xr.open_dataset(base_path / 'data' / 'nonuniform_mp_concentration.nc').mp_concentration
            mp_conc = target_density * mp_conc / mp_conc.weighted(np.cos(np.deg2rad(mp_conc.latitude))).mean()
            horizontal = HorizontalDistribution.from_dataarray(mp_conc)
    
        else:
            print(f'No horizontal distribution function exists for the option `{horizontal}`; skipping `{experiment}`')
            continue
    
        # ---------------------------------------------------------------------------------------------------------
    
        print(f'Calculating files for `{experiment}`')
        dir.mkdir(parents=True, exist_ok=True)
    
        n_levels = vertical.shape[0] if len(vertical.shape)==3 else vertical.shape[1]
    
        if len(horizontal.shape) == 3 and horizontal.shape[0] == len(time):
            horizontal = np.repeat(horizontal[:,np.newaxis,:], n_levels, axis=1)
        else:
            horizontal = np.repeat(horizontal[np.newaxis,:], n_levels, axis=0)
            horizontal = np.repeat(horizontal[np.newaxis,:], len(time), axis=0)
    
        if not len(vertical.shape) == 4 and vertical.shape[0] == len(time):
            vertical = np.repeat(vertical[np.newaxis,:], len(time), axis=0)
    
        total_distribution = vertical * horizontal
    
        for file, band in [(sw_file, 'sw'), (lw_file, 'lw')]:
            for name in ['asymmetry', 'extinction', 'absorption']:
                output_file = dir / f'easy_{name}_{band}.nc'
                if not output_file.exists():
                    data = make_3d_scattering_data(name, total_distribution, file)
                    save_netcdf(output_file, data, f'easy_{name}_{band}', '1' if name == 'asymmetry' else 'm-1')
                else:
                    print(f'  {band} {name} already present, skipping')
