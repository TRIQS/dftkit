from h5 import HDFArchive
import numpy as np


def assert_h5_close(actual_path, expected_path, atol=1e-12):
    """Compare every HDF value at the requested tolerance."""

    def compare(actual, expected, path=''):
        if hasattr(actual, 'keys'):
            assert set(actual.keys()) == set(expected.keys()), path
            for key in actual.keys():
                compare(actual[key], expected[key], f'{path}/{key}')
        elif isinstance(actual, (list, tuple)):
            assert len(actual) == len(expected), path
            for i, (left, right) in enumerate(zip(actual, expected)):
                compare(left, right, f'{path}[{i}]')
        elif isinstance(actual, (np.ndarray, np.number, int, float, complex)):
            assert np.shape(actual) == np.shape(expected), path
            dtype = np.asarray(actual).dtype
            if np.issubdtype(dtype, np.integer) or dtype.kind in 'OUS':
                np.testing.assert_array_equal(actual, expected, err_msg=path)
            else:
                assert np.isfinite(actual).all() and np.isfinite(expected).all(), path
                np.testing.assert_allclose(actual, expected, atol=atol, rtol=0,
                                           equal_nan=False, err_msg=path)
        else:
            assert actual == expected, path

    with HDFArchive(actual_path, 'r') as actual, \
            HDFArchive(expected_path, 'r') as expected:
        compare(actual, expected)
