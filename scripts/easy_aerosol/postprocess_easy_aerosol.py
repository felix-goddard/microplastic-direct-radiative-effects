import subprocess
from pathlib import Path

start_year = 1990
end_year = 2009 # inclusive
ea_dir = '/path/to/model/output'

scenarios = [
    'control',
    'clear_gamma(1-100)_exponential(0.3,2km)_uniform(420m-3)',
    'clear_gamma(1-100)_exponential(0.3,2km)_uniform(840m-3)',
    'clear_gamma(1-100)_exponential(0.3,2km)_uniform(1680m-3)',
    'clearUV_gamma(1-100)_exponential(0.3,2km)_uniform(1760m-3)',
    'black_gamma(1-100)_exponential(0.3,2km)_uniform(860m-3)',
    'mixedUV_gamma(1-100)_exponential(0.3,2km)_uniform(7230m-3)',
    'clear_powerlaw(1.52,1-100)_exponential(0.3,2km)_uniform(5710m-3)',
    'clear_gamma(1-100)_wellmixedBL_uniform(2620m-3)',
    'clear_gamma(1-100)_exponential(0.3,2km)_evangeliou(3070m-3)',
    'mixedUV_powerlaw(1.52,1-100)_wellmixedBL_evangeliou(650m-3)',
]

#==================================================================

base_path = Path(__file__).resolve().parent.parent.parent
postprocess_cmd = f'QT_QPA_PLATFORM=offscreen python {base_path}/scripts/calc_rad_eff.py'

ea_dir = Path(ea_dir)
for scenario in scenarios:
    if scenario != 'control':
        (ea_dir / 'annual' / scenario / 'data').mkdir(parents=True, exist_ok=True)
        (ea_dir / 'annual' / scenario / 'plots').mkdir(parents=True, exist_ok=True)

    for year in range(start_year, end_year+1):
        year_dir = ea_dir / 'annual' / scenario / f'{year:04}'
        year_dir.mkdir(parents=True, exist_ok=True)
        for file in (ea_dir / scenario).glob(f'{year:04}-*'):
            if not (year_dir / file.name).exists():
                (year_dir / file.name).symlink_to(file)

        if scenario != 'control':
            control_dir = ea_dir / 'annual' / 'control' / f'{year:04}'
            output_file = ea_dir / 'annual' / scenario / 'data' / f'{year:04}.nc'
            plot_file = ea_dir / 'annual' / scenario / 'data' / f'{year:04}.png'

            if not output_file.exists():
                subprocess.run(
                    f'{postprocess_cmd} "{control_dir}" "{year_dir}" "{output_file}" "{plot_file}"',
                    shell=True)
        
                print(f'Created {output_file}')
