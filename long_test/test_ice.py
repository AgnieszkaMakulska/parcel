"""
This test runs the parcel model using blk_1m_ice microphysics scheme.
It checks that the final values of ria and rc match their reference values (are consistent in each run).
"""

import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")
from parcel import parcel
from scipy.io import netcdf
import numpy as np

def test_bulk_ice():
    args = dict(
        dt=0.1,
        z_max=500,
        w=1.0,
        T_0=273,
        p_0=101300,
        RH_0 = 1,
        outfile="test_bulk_ice.nc",
        outfreq=100,
        scheme="blk_1m_ice"
    )
    parcel(**args)
    with netcdf.netcdf_file("test_bulk_ice.nc", 'r') as f:
        ria = np.array(f.variables['ria'][:])
        rc = np.array(f.variables['rc'][:])
    assert np.isclose(ria[-1], 7.843e-11, rtol=1e-3)
    assert np.isclose(rc[-1], 5.522e-4, rtol=1e-3)
