import re
import scipy
import numpy as np
import xarray as xr
from abc import ABC, abstractmethod, abstractproperty
from pathlib import Path
from . import constants
from . import size_distribution
from .grid_utils import spatial_mean, spatial_integral, EARTH_RADIUS

EARTH_SURFACE_AREA = 4 * np.pi * EARTH_RADIUS**2

base_path = Path(__file__).resolve().parent.parent
profiles_dir = base_path / 'data' / 'socrates' / 'atmospheric_profiles'


class BaseExperiment():
    def __init__(self, label):
        colour, size_dist, vertical_dist, horizontal_dist = label.split('_')

        self.size_dist = size_dist
        if (match := re.match(r'gamma\((\d+(?:\.\d+)?)-(\d+(?:\.\d+)?)\)', self.size_dist)) is not None:
            self.size_distribution = size_distribution.Revell2021GammaSizeDistribution(
                float(match.group(1))/2, float(match.group(2))/2)

        elif (match := re.match(r'powerlaw\((\d+(?:\.\d+)?),(\d+(?:\.\d+)?)-(\d+(?:\.\d+)?)\)', self.size_dist)) is not None:
            self.size_distribution = size_distribution.PowerLawSizeDistribution(
                -float(match.group(1)), float(match.group(2))/2, float(match.group(3))/2)

        else:
            raise Exception(f'Unknown size distribution `{self.size_dist}`')

        self.horizontal_dist = horizontal_dist
        if (match := re.match(r'(uniform|evangeliou)\((\d+(?:\.\d+)?)m-3\)', horizontal_dist)) is not None:
            self.horizontal_dist_type = match.group(1)
            self.surface_concentration = float(match.group(2))
        else:
            raise Exception(f'Unknown horizontal distribution `{self.horizontal_dist}`')

        self.colour = colour
        self.vertical_dist = vertical_dist
        self._mass_distribution = None

    @property
    def name(self):
        return f'{self.colour}_{self.size_dist}_{self.vertical_dist}_{self.horizontal_dist}'

    @abstractproperty
    def _lw_forcing(self):
        pass

    @abstractproperty
    def _sw_forcing(self):
        pass

    def longwave_forcing                 (self):            return self._lw_forcing

    def global_mean_longwave_forcing     (self, mask=None): return spatial_mean(self._lw_forcing, mask=mask).values

    def shortwave_forcing                (self):            return self._sw_forcing

    def global_mean_shortwave_forcing    (self, mask=None): return spatial_mean(self._sw_forcing, mask=mask).values

    def net_forcing                      (self):            return self._lw_forcing + self._sw_forcing

    def global_mean_net_forcing          (self, mask=None): return spatial_mean(self._lw_forcing + self._sw_forcing, mask=mask).values

    def longwave_efficiency              (self):            return self._lw_forcing / self.mass_distribution

    def global_mean_longwave_efficiency  (self, mask=None): return (
        spatial_mean(self._lw_forcing, mask=mask) / spatial_mean(self.mass_distribution, mask=mask)).values

    def shortwave_efficiency             (self):            return self._sw_forcing / self.mass_distribution

    def global_mean_shortwave_efficiency (self, mask=None): return (
        spatial_mean(self._sw_forcing, mask=mask) / spatial_mean(self.mass_distribution, mask=mask)).values

    def net_efficiency                   (self):            return (self._lw_forcing + self._sw_forcing) / self.mass_distribution

    def global_mean_net_efficiency       (self, mask=None): return (
        spatial_mean(self._lw_forcing + self._sw_forcing, mask=mask) / spatial_mean(self.mass_distribution, mask=mask)).values

    @property
    def mass_distribution(self):
        """Get the distribution of plastic mass for this experiment."""

        if self._mass_distribution is not None:
            return self._mass_distribution

        # load the profiles of temperature and water vapour
        layer_t = (xr.open_dataset(profiles_dir / 'temperature_layers.nc').t.transpose('plev','lat','lon'))
        center_t = (xr.open_dataset(profiles_dir / 'temperature_centers.nc').t.transpose('plev','lat','lon'))
        center_q = (xr.open_dataset(profiles_dir / 'water_vapour.nc').q.transpose('plev','lat','lon'))
    
        # calculate the atmospheric layer thicknesses
        virtual_temperature = center_t * (1 + center_q * (1 / constants.water_vapour_mass_ratio - 1))
        log_pressure_diff = np.log(layer_t.plev).diff('plev').assign_coords(plev=layer_t.plev-layer_t.plev.diff('plev')/2)
        thickness = (constants.dry_gas_constant / constants.gravity) * virtual_temperature * log_pressure_diff
    
        # calculate the average particle mass
        expected_volume = self.size_distribution.nth_moment(3) * (4/3) * np.pi
        particle_mass = constants.plastic_density * expected_volume * 1e-18 # convert um^3 to m^3 since the size distributions work in um

        self._mass_distribution = (particle_mass * self.number_concentration * thickness).sum('plev')
        return self._mass_distribution

    @property
    def number_concentration(self):
        filename = f"mps_{self.vertical_dist}_{self.horizontal_dist_type}(1m-3).nc"
        mp_concentration = xr.open_dataset(profiles_dir / filename).mp_conc
        return mp_concentration * self.surface_concentration

    @property
    def total_plastic_burden(self):
        return spatial_integral(
            self.mass_distribution.transpose('lat','lon'),
            longitude_name='lon', latitude_name='lat').values

    def concentration_for_forcing_magnitude(self, forcing, mask=None):
        efficiency = self.global_mean_net_efficiency(mask=mask)
        surface_concentration = self.surface_concentration
        burden_per_conc = self.total_plastic_burden / surface_concentration
        return abs(forcing / efficiency) / (burden_per_conc / EARTH_SURFACE_AREA)


#====================================================================================================================================


class EasyAerosolExperiment(BaseExperiment):
    def __init__(self, label):
        super().__init__(label)

        folder_name = f'{self.colour}_{self.size_dist}_{self.vertical_dist}_{self.horizontal_dist}'
        
        path = base_path / 'data' / 'easy_aerosol' / 'output' / folder_name
        self._data = xr.concat([
            xr.open_dataset(f, engine='h5netcdf').expand_dims(dim={'year': [int(f.stem)]}, axis=0)
            for f in path.glob('*.nc')
        ], dim='year').sortby('year')

        self._mass_distribution = None

    @property
    def _lw_forcing(self):
        return self._data.rlut_cs2

    @property
    def _sw_forcing(self):
        return self._data.rsut_cs2

    @staticmethod
    def tdistribution_95ci(x):
        n = len(x)
        mean = np.mean(x)
        err = scipy.stats.t.isf((100-90)/100., n - 1, 0, np.std(x)/np.sqrt(n))
        return np.array([mean-err, mean, mean+err])

    def longwave_forcing_95ci(self, mask=None):
        ts = self.global_mean_longwave_forcing(mask=mask)
        return self.tdistribution_95ci(ts)

    def longwave_efficiency_95ci(self, mask=None):
        ts = self.global_mean_longwave_efficiency(mask=mask)
        return self.tdistribution_95ci(ts)

    def shortwave_forcing_95ci(self, mask=None):
        ts = self.global_mean_shortwave_forcing(mask=mask)
        return self.tdistribution_95ci(ts)

    def shortwave_efficiency_95ci(self, mask=None):
        ts = self.global_mean_shortwave_efficiency(mask=mask)
        return self.tdistribution_95ci(ts)

    def net_forcing_95ci(self, mask=None):
        ts = self.global_mean_net_forcing(mask=mask)
        return self.tdistribution_95ci(ts)

    def net_efficiency_95ci(self, mask=None):
        ts = self.global_mean_net_efficiency(mask=mask)
        return self.tdistribution_95ci(ts)


#====================================================================================================================================


class SocratesExperiment(BaseExperiment):
    def __init__(self, label):
        super().__init__(label)

        expt = f'{self.colour}_{self.size_dist}_{self.vertical_dist}_{self.horizontal_dist}'

        output_dir = base_path / 'data' / 'socrates' / 'output'
        assert Path(output_dir / f'{expt}.nc').exists(), f'No data exists for experiment `{expt}`'
        
        data = xr.open_dataset(output_dir / f'{expt}.nc')
        control = xr.open_dataset(output_dir / 'control.nc')

        self.__lw_forcing = (data.lw - control.lw).isel(plev=0)
        self.__sw_forcing = (data.sw - control.sw).isel(plev=0)

        self._mass_distribution = None

    @property
    def _lw_forcing(self):
        return self.__lw_forcing

    @property
    def _sw_forcing(self):
        return self.__sw_forcing



