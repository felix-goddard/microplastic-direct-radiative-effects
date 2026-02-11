# Revised direct radiative forcing of airborne microplastics suggests warming

Felix W. Goddard<sup>1,∗</sup>, Stefania Glukhova<sup>2</sup>, Eric C. Le Ru<sup>2</sup>, Nikolaos Evangeliou<sup>3</sup>, Cameron McErlich<sup>1</sup>, Catherine Hardacre<sup>1</sup>, Dave Frame<sup>1</sup>, Peter Kuma<sup>4</sup>, and Laura E. Revell<sup>1,∗</sup>

<sup>1</sup> School of Physical and Chemical Sciences, University of Canterbury, Christchurch, New Zealand<br>
<sup>2</sup> The MacDiarmid Institute for Advanced Materials and Nanotechnology, School of Chemical and Physical Sciences, Victoria University of Wellington, Wellington, New Zealand<br>
<sup>3</sup> Stiftelsen NILU (formerly NILU - Norwegian Institute for Air Research), Department for Atmospheric & Climate Research (ATMOS), 2007 Kjeller, Norway<br>
<sup>4</sup> Rossby Centre, Swedish Meteorological and Hydrological Institute, Norrköping, Sweden<br>
<sup>*</sup> Correspondence to: felix.goddard@pg.canterbury.ac.nz and laura.revell@canterbury.ac.nz

This repository contains software and data accompanying the manuscript by Goddard et al. (2026) "Revised direct radiative forcing of airborne microplastics suggests warming".

# Requirements

The software should work as intended on any Linux or macOS operating system; it has been tested on Rocky Linux version 9.4 with Python 3.11.9. Additional packages required are available on the Python Package Index; they are:
- cartopy 0.25.0
- cdsapi 0.7.7
- cftime 1.6.4
- ds-format 1.1.0
- h5netcdf 1.7.3 (optional; netCDF4 will work in its place, but will spew warnings due to the use of parentheses in filenames)
- jupyter 1.1.1
- matplotlib 3.10.7
- miepython 2.5.4
- netcdf4 1.7.3
- numpy 2.3.4
- pandas 2.3.3
- scicomap 1.0.1
- scipy 1.16.3
- setuptools 80.9.0 (specific version required by scicomap)
- xarray 2025.10.1

A conda environment containing these packages can be created using the provided `environment.yml` file:
```shell
$ conda env create -f environment.yml
```

# Figures

All figure generation code is contained in the notebook `Figures.ipynb`; this notebook operates only on data included in this repository, so should work out-of-the-box. The whole notebook can be run in a few minutes.