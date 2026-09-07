################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Aichhorn, L. Pourovskii, V. Vildosola
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

"""
Test the high-symmetry k-path label block that dmftproj appends to case.outband.
A real case.outband is far too large to ship as a fixture, so this writes a
synthetic tail in the exact Fortran (2i6,a) format of outband.f:277, including
the trailing blank the `a` descriptor emits for the declared length of kname.
"""

import os
import tempfile

import numpy as np

from triqs_dftkit.wien2k.converter import _read_kpath_labels

N_K = 501


def line(counter, pos, label):
    return '%6d%6d%s' % (counter, pos, label.ljust(10))


def read(lines, n_k=N_K):
    fd, path = tempfile.mkstemp(suffix='.outband')
    try:
        with os.fdopen(fd, 'w') as f:
            f.write('\n'.join(lines) + '\n')
        return _read_kpath_labels(path, n_k)
    finally:
        os.remove(path)


# One block exercising the awkward cases at once: '\xG' does not start with a
# letter (XCrySDen writes Gamma that way), X|Y puts two labels on the same
# k-point, the projector line in front of the block itself parses as (2i6,a),
# the padding pushes the block past the tail that is read back, and the file
# ends on blank lines.
padding = ['%20.14f%20.14f' % (0.5, 0.25)] * 2000
block = [line(1, 1, '\\xG'), line(2, 122, 'X'), line(3, 122, 'Y'),
         line(4, 268, 'L'), line(5, 501, 'K')]
labels, idx = read(padding + [line(3, 7, '0.52')] + block + ['', '   '])

assert labels == ['\\xG', 'X', 'Y', 'L', 'K'], labels
np.testing.assert_array_equal(idx, [0, 121, 121, 267, 500])

# A malformed line in the middle of the block truncates the backward walk, so
# the block no longer starts at counter 1 and is rejected rather than silently
# returning only its tail.
broken = list(block)
broken[2] = 'this line is not a label line'
assert read(padding + broken) == (None, None)

# A position past the end of the band path is rejected as well.
assert read(padding + block, n_k=400) == (None, None)
