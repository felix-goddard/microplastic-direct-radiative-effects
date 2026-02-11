This directory contains scripts for generating the input data and running SOCRATES.

The basic workflow using these scripts is as follows:
- `create_profiles.py` generates netCDF files describing the atmospheric structure that SOCRATES will simulate;
- `create_spectra.py` reads the optical property data used by EasyAerosol (see `scripts/easy_aerosol`) and formats it in a manner SOCRATES can ingest;
- `run_socrates` runs SOCRATES.

# `create_profiles.py`

This script reads in data that has been downloaded and regridded from ERA5 and formats it in the manner expected by SOCRATES. This can simply be run:
```shell
$ python create_profiles.py
```

The ERA5 input data to this script is the files `data/socrates/atmos_climatology.nc` and `data/socrates/surface_climatology.nc`. These were downloaded using the [cdsapi package](https://github.com/ecmwf/cdsapi); see `get_era5.py`. The downloaded data was postprocessed by time-averaging and regridding to the UM grid using CDO and NCO, e.g. for the atmospheric data:
```bash
# rename the time dimension from "date" to just "time"
ncrename -d date,time -v date,time downloaded_file.nc tmp.nc
# # select variables, take a time average, regrid to the UM grid
cdo -remapcon,umgridfile -timmean -selname,t,q,o3 tmp.nc atmos_climatology.nc
```

The file describing the UM grid `umgridfile` can be found in `data/socrates/`.

Each profile for microplastics is only produced at a global-average concentration of 1 MP m<sup>-3</sup>. The script that calls SOCRATES itself handles scaling this profile to the specified concentration.

# `create_spectra.py`

This script reads the optical property data produced for EasyAerosol and formats it in the manner expected by SOCRATES. This requires the spectra files provided with the SOCRATES code: `sp_sw_ga7` in the shortwave and `sp_lw_ga7` in the longwave. Line 3 should be modified to point to the path to the base of the SOCRATES directory tree.

Once the script is set up to point to the SOCRATES code, it is configured to run on all the files of optical properties it finds in `data/optical_properties`, so can simply be run:
```shell
$ python create_spectra.py
```
This will create spectra files in `data/socrates/optical_properties`.

As SOCRATES is designed to handle a pre-determined set of aerosol species, I have set this up to co-opt the "soot" aerosol species; the atmospheric profiles produced by the `create_profiles.py` script reflect this too.

# `run_socrates`

This script handles the setup and running of SOCRATES itself. Line 3 defines a variable `SOCRATES_PATH` which should be filled with the path to the base of the SOCRATES directory tree.

Lines 9--12 of the script define the parameters of the simulated microplastics. These correspond to the four terms of the EasyAerosol experiment names, described in `scripts/easy_aerosol/README.md`, and allow for the same parameter settings. If a combination of vertical profile (set by `mp_vertical`) and horizontal profile (set by `mp_horizontal`) is not valid, the script will default to running with no microplastics (i.e. the control configuration).

As loading data and calculating radiative transfer for the entire globe at once requires too much memory, the script splits the globe into a series of tiles and simulates each individually before stitching the results back together. This makes use of CDO and NCO.

This will output netCDF files in `data/socrates/output`. These files have variables `sw` and `lw` corresponding to the shortwave and longwave radiative fluxes respectively.
