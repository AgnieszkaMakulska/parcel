"""
This test runs the parcel model using three different microphysics schemes:
- lgrngn (Lagrangian particle-based)
- blk_1m (bulk warm)
- blk_1m_ice (bulk ice)

It compares the evolution of rv, th_d, and total condensed water in different schemes. 
"""

import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")

import numpy as np
from parcel import parcel
from scipy.io import netcdf
import matplotlib.pyplot as plt

def run_scheme(scheme, outfile):
    args = dict(
        dt=0.1,
        z_max=800,
        w=1.0,
        T_0=300,
        r_0=0.022,
        outfile=outfile,
        outfreq=50,
        scheme=scheme,
        out_bin='{"radius": {"rght": 1, "moms": [3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 1e-15}}'
    )
    parcel(**args)
    with netcdf.netcdf_file(outfile, 'r') as f:
        rv = np.array(f.variables['r_v'][:])
        th_d = np.array(f.variables['th_d'][:])
        z = np.array(f.variables['z'][:])
        if scheme.startswith("blk"):
            r_tot = np.array(f.variables['rc'][:]) + np.array(f.variables['rr'][:])
        else: 
            moment_3 = np.array(f.variables['radius_m3'][:])
            r_tot = moment_3 *4/3 * np.pi * 997 #multiply by density of water
    return rv, th_d, z, r_tot

def test_compare_schemes():
    schemes = ["lgrngn", "blk_1m", "blk_1m_ice"]
    results = {}
    for scheme in schemes:
        rv, th_d, z, r_tot = run_scheme(scheme, f"test_{scheme}.nc")
        results[scheme] = (rv, th_d, z, r_tot)
    # Compare final values
    rv_vals = [results[s][0][-1] for s in schemes]
    th_d_vals = [results[s][1][-1] for s in schemes]
    r_tot_vals = [results[s][2][-1] for s in schemes]

    fig, ax = plt.subplots(1,3, figsize=(12, 6))
    for scheme in schemes:
        rv, th_d, z, r_tot = results[scheme]
        ax[0].plot(rv, z, label=f"{scheme}")
        ax[1].plot(th_d, z, label=f"{scheme}")
        ax[2].plot(r_tot, z, label=f"{scheme}")
    ax[0].set_ylabel("Height [m]")
    ax[0].set_xlabel("Water vapor mixing ratio")
    ax[1].set_xlabel("Dry potential temperature")
    ax[2].set_xlabel("Total condensed water mixing ratio")
    ax[0].legend()
    ax[1].legend()
    ax[2].legend()
    plt.tight_layout()
    plt.savefig("plots/outputs/scheme_comparison.svg")
    
    # Check closeness
    for i in range(1, len(schemes)):
      assert np.isclose(rv_vals[0], rv_vals[i], rtol=5e-4), f"r_v differs: {rv_vals[0]} vs {rv_vals[i]}"
      assert np.isclose(th_d_vals[0], th_d_vals[i], rtol=5e-4), f"th_d differs: {th_d_vals[0]} vs {th_d_vals[i]}"
