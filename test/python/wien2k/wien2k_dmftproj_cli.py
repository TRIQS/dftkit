################################################################################
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
# (GPL-3.0-or-later)
################################################################################

import gzip
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


here = Path(__file__).resolve().parent
wrapper = Path(sys.argv[1]).resolve()
with tempfile.TemporaryDirectory() as tmp:
    case_dir = Path(tmp) / 'CaOs2_noso_partial'
    case_dir.mkdir()
    case = case_dir.name
    shutil.copy(here / f'{case}.indmftpr', case_dir / f'{case}.indmftpr')
    for ext in ('struct', 'dmftsym'):
        shutil.copy(here / f'CaOs2.{ext}', case_dir / f'{case}.{ext}')
    for spin in ('up', 'dn'):
        source = here / f'CaOs2_noso.almblm{spin}.gz'
        with gzip.open(source, 'rb') as compressed:
            (case_dir / f'{case}.almblm{spin}').write_bytes(compressed.read())
    run = subprocess.run([sys.executable, str(wrapper), '-sp'], cwd=case_dir,
                         env=os.environ, stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT, text=True)
    assert run.returncode == 0, run.stdout
    for suffix in ('ctqmcout', 'symqmc', 'sympar', 'parproj',
                   'oubwinup', 'oubwindn'):
        output = case_dir / f'{case}.{suffix}'
        assert output.is_file() and output.stat().st_size > 0, output

print('wien2k_dmftproj_cli: ok')
