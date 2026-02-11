import os
import sys
import miepython
import numpy as np
import pandas as pd
from pathlib import Path

# Add the lib folder to the path so we can import size_distribution
base_path = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(base_path))

import lib.size_distribution as size_distribution


refractive_indices_dir = base_path / 'data' / 'refractive_indices'
spectra_dir = base_path / 'data' / 'spectra'


#==============================================================================================================


def calculate_mie_properties(refractive_indices, size_distribution, n_radii=100):
    '''Calculates size-averaged Mie properties (extinction and scattering cross-sections
    and asymmetry parameter.
    
    refractive_indices should be a Pandas dataframe with columns for wavelength (in microns)
      and the real (n) and imaginary (k) parts of the refractive index.
      
    size_distribution should be a SizeDistribution (defined in lib/size_distribution.py)
    
    n_radii defines the number of bins of particle radius over which we average the spectra.'''
    
    wavelength = refractive_indices.wavelength
    refractive_index = refractive_indices.n - 1j * refractive_indices.k

    qext = np.zeros_like(wavelength)
    qsca = np.zeros_like(wavelength)
    g = np.zeros_like(wavelength)
    total_weight = 0

    particle_radius_edges = np.logspace(
        np.log10(size_distribution.r_min), np.log10(size_distribution.r_max), n_radii+1)

    d_radius = np.diff(particle_radius_edges)
    particle_radius_centers = (particle_radius_edges[1:] + particle_radius_edges[:-1]) / 2

    for radius, dr in zip(particle_radius_centers, d_radius):
        area = np.pi * radius**2
        size_parameter = np.pi * 2*radius / wavelength
        weight = size_distribution.pdf(radius) * dr
        ext, sca, _, asy = miepython.mie(refractive_index, size_parameter)
        qext += ext * area * weight
        qsca += sca * area * weight
        g += asy * weight
        total_weight += weight

    return qext/total_weight, qsca/total_weight, g/total_weight


def load_refractive_indices(file):
    refractive_index = pd.read_csv(refractive_indices_dir / f'{file}.csv')
    refractive_index['wavelength'] = refractive_index.wavelength.map(lambda x: x/1000) # Convert wavelength from nm to μm
    return refractive_index


def create_spectra(refractive_index, size_dist):
    ext, sca, asy = calculate_mie_properties(refractive_index, size_dist)
    return pd.DataFrame(
        data={'wavelength [μm]': refractive_index.wavelength,
              'Cabs [μm^2]': ext-sca,
              'Csca [μm^2]': sca,
              'g [1]': asy}
    ).set_index('wavelength [μm]')


#==============================================================================================================


if __name__ == '__main__':
    for colour in ['clear', 'clearUV', 'mixedUV', 'black']:
        refractive_index = load_refractive_indices(colour)
    
        size_dist = size_distribution.Revell2021GammaSizeDistribution(.5, 50)
        create_spectra(refractive_index, size_dist).to_csv(spectra_dir / f'{colour}_gamma(1-100).csv')
    
        size_dist = size_distribution.PowerLawSizeDistribution(-1.52, .5, 50)
        create_spectra(refractive_index, size_dist).to_csv(spectra_dir / f'{colour}_powerlaw(1.52,1-100).csv')

    for colour in ['blue365', 'indigo2008', 'pink4815', 'pink4890']:
        refractive_index = load_refractive_indices(colour)
    
        size_dist = size_distribution.Revell2021GammaSizeDistribution(.5, 50)
        create_spectra(refractive_index, size_dist).to_csv(spectra_dir / f'{colour}_gamma(1-100).csv')
