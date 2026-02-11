from pathlib import Path

socrates_path = '/path/to/socrates'

#=====================================================

base_path = Path(__file__).resolve().parent.parent.parent

src_dir = base_path / 'data' / 'optical_properties'
out_dir = base_path / 'data' / 'socrates' / 'optical_properties'

out_dir.mkdir(parents=True, exist_ok=True)

sw_start_index = 387 # in the sp_sw_ga7 file, this is the line index where we need to fill from
lw_start_index = 560 # likewise for sp_lw_ga7

spectra_dir = Path(socrates_path) / 'data' / 'spectra' / 'ga7'

sw_template = None
with open(spectra_dir / 'sp_sw_ga7') as f:
    sw_template = f.readlines()

lw_template = None
with open(spectra_dir / 'sp_lw_ga7') as f:
    lw_template = f.readlines()

files = [x for x in src_dir.iterdir() if x.is_file()]
for file in files:
    output_file = out_dir / file.name
    if output_file.exists():
        continue

    with open(file) as f:
        data_lines = [l for l in f.readlines()[13:] if l.strip()]
        file_lines = (sw_template if 'sw' in file.name else lw_template if 'lw' in file.name else None).copy()
        start_idx = sw_start_index if 'sw' in file.name else lw_start_index if 'lw' in file.name else None
        for i, line in enumerate(data_lines):
            file_lines[start_idx+i] = line
        with open(output_file, 'w') as g:
            g.writelines(file_lines)

    print('Created', file.name)