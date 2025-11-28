"""
This test runs the parcel model using lagrangian ice microphysics.
It checks that the final results match their reference values (are consistent in each run).
"""

import sys, os, subprocess
sys.path.insert(0, "../../")
sys.path.insert(0, "./")
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from parcel import parcel
from scipy.io import netcdf


def test_ice_SD():    
    ref_dir = "long_test/refdata/"
    outfile = "test.nc"
    for rng_seed in [1234, 5678, 91011]:
       for (ref_file, time_dep) in (ref_dir+"ice_ref_sing.nc", False), (ref_dir+"ice_ref_timedep.nc", True):
              parcel(dt=1.,w=1.,sd_conc=100,
                     z_max = 5000.0,
                     T_0 = 263.0,
                     RH_0 = 1.,
                     scheme = "lgrngn",
                     ice_switch=True,
                     ice_nucl=True,
                     time_dep_ice_nucl=time_dep,
                     aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "rd_insol": 0.5e-6, "mean_r": [0.02e-6], "gstdev": [1.4], "n_tot": [60.0e6]}}', 
                     outfreq = 5000, 
                     out_bin= '{"liq": {"rght": 1, "moms": [0,3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 5e-20}}',
                     outfile=outfile,
                     rng_seed=int(rng_seed))
              fnc = netcdf.netcdf_file(outfile)
              fnc_ref = netcdf.netcdf_file(ref_file)
              for variable in ["r_v", "th_d", "ice_mass"]:
                     assert np.isclose(fnc.variables[variable][:][-1], fnc_ref.variables[variable][:][-1], rtol=1)
              fnc.close()
              fnc_ref.close()
              subprocess.call(["rm", outfile])