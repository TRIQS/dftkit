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

from triqs_dftkit.wien2k import Driver


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

    # The CSC workflow uses the Python API directly by default. A dummy WIENROOT
    # is enough here because this focused call must not launch the x script.
    wienroot = case_dir / 'wienroot'
    wienroot.mkdir()
    (wienroot / 'x').write_text('must not be executed\n')
    previous_cwd = Path.cwd()
    try:
        os.chdir(case_dir)
        Driver(case, wienroot=str(wienroot))._run_dmftproj()
    finally:
        os.chdir(previous_cwd)

print('wien2k_dmftproj_cli: ok')
