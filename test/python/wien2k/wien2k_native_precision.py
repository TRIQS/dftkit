"""Check native precision against unitarity, independently of golden outputs."""

import gzip
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

import numpy as np

here = Path(__file__).resolve().parent
with tempfile.TemporaryDirectory() as tmp:
    case_dir = Path(tmp) / 'CaOs2'
    case_dir.mkdir()
    for suffix in ('struct', 'dmftsym', 'indmftpr'):
        shutil.copyfile(here / f'CaOs2.{suffix}', case_dir / f'CaOs2.{suffix}')
    for spin in ('up', 'dn'):
        with gzip.open(here / f'CaOs2.almblm{spin}.gz', 'rb') as src:
            (case_dir / f'CaOs2.almblm{spin}').write_bytes(src.read())
    templates = Path(tmp) / 'SRC_templates'
    templates.mkdir()
    # Wien2k supplies the d template externally. Use full precision here to
    # distinguish constructor rounding from truncation in the input template.
    basis = (here / 'case.cf_d_eg_t2g').read_text()
    basis = basis.replace('.70710678', '.7071067811865475244')
    (templates / 'case.cf_d_eg_t2g').write_text(basis)
    run = subprocess.run([sys.argv[1], '-so', '-sp'], cwd=case_dir,
                         env=dict(os.environ, WIENROOT=tmp),
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         text=True)
    assert run.returncode == 0, run.stdout
    tokens = (case_dir / 'CaOs2.symqmc').read_text().split()
    header_size = 2 + 16 * 3 + 16
    blocks = np.array(tokens[header_size:], dtype=float).reshape(16, 2, 2, 10, 10)
    matrices = blocks[:, :, 0] + 1j * blocks[:, :, 1]
    for operation in matrices:
        for matrix in operation:
            np.testing.assert_allclose(matrix @ matrix.conj().T, np.eye(10),
                                       atol=1e-12, rtol=0)

print('wien2k_native_precision: ok')
