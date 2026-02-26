This directory contains scripts for generating the input files for EasyAerosol.

The basic workflow using these scripts is as follows:
- `generate_spectra.py` is used to calculate Mie scattering and absorption spectra from the measured refractive indices (in `data/refractive_indices`);
- `prepare_cscatter.py` reads the CSV files containing the Mie spectra and formats them in the manner expected by SOCRATES;
- `run_cscatter_average` invokes the SOCRATES `Cscatter_average` script to average the Mie spectra  across the model bands, weighted by the solar spectrum (Lean & DeLand, 2012).
- `create_easy_aerosol_input.py` generates the EasyAerosol input files.

# `generate_spectra.py`

This script takes in the CSV files containing the measured plastic complex refractive indices (in `data/refractive_indices`) and calculates, for a given size distribution (see `lib/size_distribution.py`), the Mie scattering and absorption spectra and asymmetry factor, averaged across the size distribution. It can simply be run:
```shell
$ python generate_spectra.py
```

Lines 79 onwards define the spectra that will be generated. The general formula is
```python
# `colour` is one of the files in data/refractive_indices (without the .csv extension)
refractive_index = load_refractive_indices(colour)

# an instance of one of the child classes of size_distribution.SizeDistribution, e.g.
size_distribution = size_distribution.Revell2021GammaSizeDistribution(.5, 50)

spectra = create_spectra(refractive_index, size_distribution)
```
The returned object `spectra` is a pandas DataFrame, which the script then saves to a CSV file in `data/spectra`.

# `prepare_cscatter.py`

This script reads the CSV files produced by `generate_spectra.py` and processes them into a format expected by SOCRATES in the following step. It is configured to work out-of-the box by reading the CSV files in the `data/spectra` and producing corresponding files lacking the `.csv` extension that SOCRATES can ingest. It can thus simply be run:
```shell
$ python prepare_cscatter.py
```

# `run_cscatter_average`

This script reads the Mie spectra produced in the previous steps and averages them across the model bands, saving the output in `data/optical_properties`. It must be run with an argument for the name of the file to be processed, e.g.
```shell
$ ./run_cscatter_average 'clear_gamma(1-100)'
```

This script invokes the SOCRATES script `Cscatter_average`, which must be available. Line 3 of the script defines a variable `SOCRATES_PATH` which should be filled with the path to the base of the SOCRATES directory tree (the script will print a warning if this isn't done).

# `create_easy_aerosol_input.py`

This script generates the shortwave and longwave absorption, extinction, and asymmetry files that drive the EasyAerosol code.

The files it generates are specified by the `experiments` variable defined on line 12. Simply running the code:
```shell
$ python create_easy_aerosol_input.py
```
will create a folder in `data/easy_aerosol/forcing` corresponding to each element of the `experiments` array and fill that folder with the EasyAerosol input files.

Each experiment is described by a name consisting of four terms separated by underscores, e.g. `clear_gamma(1-100)_exponential(0.3,2km)_uniform(1680m-3)`. The four terms are, in order: microplastic colour, particle size distribution, vertical spatial distribution, and horizontal spatial distribution.

The variable `opt_prop_dir` on line 27 specifies the path to the directory containing the files of optical properties (`data/optical_properties` by default); each file encodes both a colour and size distribution, hence the allowed values of colour and size distribution are determined by the files found therein.

The variable `levels_file` on line 38 specifies the path to a namelist file that defines the vertical levels of the HadGEM3-GA7.1. The level set we use is given (in namelist format) as L85(50$_{\rm t}$,35$_{\rm s}$)$_{85}$ in Section 3.2 of [the supplement to Walters et al. (2017)](https://gmd.copernicus.org/articles/10/1487/2017/gmd-10-1487-2017-supplement.pdf).

Two vertical distributions are defined (lines 214--222):
- the "exponential" distribution, which is parameterised by an exponential base and a cutoff altitude:
  `exponential(x,ykm)`, where `x` and `y` are numbers, yields the vertical distribution function $x^{z/10\text{ km}}$ for $z\leq y$ and $0$ for $z>y$.
- the "well-mixed boundary layer" distribution: `wellmixedBL`.

Two horizontal distributions are defined (lines 226--237):
- the "uniform" distribution, which is parameterised by the surface concentration: `uniform(xm-3)`, where `x` is a number, yields a
  constant surface concentration of `x` particles per cubic meter.
- the "evangeliou" distribution, based on the time-average surface microplastic concentration of a 10-year simulation of the UKESM1.1 with microplastics added as an aerosol species (McErlich et al. 2025) driven by the emissions inventory of Tichý et al. (2025). This distribution is parameterised by the global-average surface concentration: `evangeliou(xm-3)`, where `x` is a number, yields the non-uniform surface distribution scaled to a global average surface microplastic concentration of `x` particles per cubic meter.

# `postprocess_easy_aerosol.py`

As we have configured HadGEM3-GA7.1, the model outputs netCDF files of daily averages with the naming scheme `YYYY-MM-DD.nc`. This script reads those files and produces annual-mean netCDF files for the radiative flux variables.

The script expects the following directory structure: line 6 specifies a base directory containing folders corresponding to each experiment; the variable `scenarios` on line 8 gives the names of each of these folders. Each of those folders should contain the daily average netCDF files output by the model. The script will create a folder `annual` and then subfolders corresponding to each of the experiments, and fill those with the postprocessed data:
```
path / to / model / data
 ├─ expt_1
 │   ├─ 1990-01-01.nc
 │   ├─ 1990-01-02.nc
 │   └─ ...
 ├─ expt_2
 ├─ expt_3
 ├─ ...
 └─ annual
     ├─ expt_1
     ├─ expt_2
     └─ ...
```

Once the script has been modified to point to the correct directories, it can simply be run:
```shell
$ python postprocess_easy_aerosol.py
```

This script invokes the script `calc_rad_eff.py` which reads the model output and calculates the averages, and which is re-used without modification from [the software accompanying Revell et al. (2021)](https://github.com/peterkuma/microplastics2021/blob/master/bin/calc_rad_eff).

# References

Lean, J. L., & DeLand, M. T. (2012). How Does the Sun’s Spectrum Vary? https://doi.org/10.1175/JCLI-D-11-00571.1

McErlich, C., Goddard, F., Aves, A., Hardacre, C., Evangeliou, N., Hewitt, A. J., & Revell, L. E. (2025). Description and evaluation of airborne microplastics in the United Kingdom Earth System Model (UKESM1.1) using GLOMAP-mode. Geoscientific Model Development, 18(22), 8827–8854. https://doi.org/10.5194/gmd-18-8827-2025

Revell, L. E., Kuma, P., Le Ru, E. C., Somerville, W. R. C., & Gaw, S. (2021). Direct radiative effects of airborne microplastics. Nature, 598(7881), 462–467. https://doi.org/10.1038/s41586-021-03864-x

Tichý, O., Košík, V., Šmídl, V., & Evangeliou, N. (2025, March 18). Atmospheric microplastics emissions estimation and uncertainty quantification using Gibbs sampler. https://doi.org/10.5194/egusphere-egu25-11924

Walters, D., Boutle, I., Brooks, M., Melvin, T., Stratton, R., Vosper, S., Wells, H., Williams, K., Wood, N., Allen, T., Bushell, A., Copsey, D., Earnshaw, P., Edwards, J., Gross, M., Hardiman, S., Harris, C., Heming, J., Klingaman, N., … Xavier, P. (2017). The Met Office Unified Model Global Atmosphere 6.0/6.1 and JULES Global Land 6.0/6.1 configurations. Geoscientific Model Development, 10(4), 1487–1520. https://doi.org/10.5194/gmd-10-1487-2017

