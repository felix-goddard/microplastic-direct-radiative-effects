import sys
import pandas as pd
from pathlib import Path

base_path = Path(__file__).resolve().parent.parent.parent
spectra_dir = base_path / 'data' / 'spectra'

header = '''
*FILE TYPE =     6

                Single scattering parameters.


Scattering algorithm =     1   (Full Mie scattering calculation         )
Scattering type      =     0   (Unassigned          )

Number density   =
Effective radius =
Volume fraction  =

    Wavelength (m)      Absorption (m-1)    Scattering (m-1)     Moments of pf'''

for file in spectra_dir.glob('*.csv'):
    if not (spectra_dir / file.stem).exists():
        # Scaling factors for each column in case its not in SI units;
        # by default we will assume the columns are in µm, µm^2 etc.
        factors = [1e-6, 1e-12, 1e-12, 1] # µm, µm^2, µm^2, 1 
        # factors = [1, 1, 1, 1] # SI units

        df = pd.read_csv(file)
        lam  = df[df.columns[0]].values * factors[0]
        cabs = df[df.columns[1]].values * factors[1]
        csca = df[df.columns[2]].values * factors[2]
        g    = df[df.columns[3]].values * factors[3]

        with open(spectra_dir / file.stem, 'w') as f:
            print(header, file=f)
            for i in range(len(lam)):
                print('%20.9E%20.9E%20.9E%20.9E' % (lam[i], cabs[i], csca[i], g[i]), file=f)
            print(file=f)

            print(f'Created {file.stem}')
