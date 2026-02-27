"""
Checking if mixing between substeps is important for the WBF process
"""

import sys, os
sys.path.insert(0, "../")
sys.path.insert(0, "./")

import numpy as np
from parcel import parcel
from scipy.io import netcdf
import matplotlib.pyplot as plt
from functions import rh_to_rh_i

sstp = 10

def run_scheme(outfile, mixing, *, sstp_cond=sstp):
    args = dict(
        p_0=100000,
        RH_0=0.9,
        T_0=275,
        aerosol = None,
        sd_conc=100,
        dt=1,
        z_max=1000,
        w = 1.,
        outfile=outfile,
        outfreq=10,
        scheme="lgrngn",
        out_bin='{"liq": {"rght": 1, "moms": [0,1,3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0.5e-10},' \
                '"ice": {"rght": 1, "moms": [0,1,3], "drwt": "ice_a", "nbin": 1, "lnli": "lin", "left": 0.5e-10}}',
        sstp_cond=sstp_cond,
        adaptive_sstp_cond=False,
        sstp_cond_mix   = mixing,
        exact_sstp_cond = True,
        aerosol_independent_of_rhod=True, 
        backend="OpenMP",
        ice_switch = True,
        ice_nucl = True,
        time_dep_ice_nucl = False,
        depo = True
    )

    if hasattr(run_scheme, "aerosol"):
        args["aerosol"] = run_scheme.aerosol  

    parcel(**args)

    with netcdf.netcdf_file(outfile, 'r') as f:
        z = np.array(f.variables['z'][:])
        RH = np.array(f.variables['RH'][:])
        T = np.array(f.variables['T'][:])
        ice_mix_ratio = np.array(f.variables['ice_mix_ratio'][:])
        liq_mix_ratio = np.array(f.variables['liq_m3'][:]) *4/3 * np.pi * 997 #multiply by density of water
        liq_conc = np.array(f.variables['liq_m0'][:])  # 1/kg
        ice_conc = np.array(f.variables['ice_m0'][:])  # 1/kg
        liq_r = np.where(liq_conc > 0, np.array(f.variables['liq_m1'][:]) / np.array(f.variables['liq_m0'][:]), 0)
        ice_r = np.where(ice_conc > 0, np.array(f.variables['ice_m1'][:]) / np.array(f.variables['ice_m0'][:]), 0)
    return RH, T, z, ice_mix_ratio, liq_mix_ratio, liq_conc, ice_conc, liq_r, ice_r



def make_figure(aerosol):
    run_scheme.aerosol = aerosol
    fig, ax = plt.subplots(1, 5, figsize=(15.0, 8.0), sharey=True, squeeze=False)

    for mixing in [True, False]:

        outfile = f"test_sstp_WBF_mix_{mixing}.nc"
        RH, T, z, ice_mix_ratio, liq_mix_ratio, liq_conc, ice_conc, liq_r, ice_r = run_scheme(outfile, mixing)

        (ice_c, liq_c, ice_l, liq_l) = ("skyblue", "coral", "ice mixing", "liquid mixing") if mixing else ("steelblue", "sienna", "ice non mixing", "liquid non mixing")

        ax[0,0].plot(ice_mix_ratio * 1e3, z, label="ice", color=ice_c)
        ax[0,0].plot(liq_mix_ratio * 1e3, z, label="liquid", color=liq_c)
        ax[0,1].plot(ice_conc / 1e6, z, label="ice", color=ice_c)
        ax[0,1].plot(liq_conc / 1e6, z, label="liquid", color=liq_c)
        ax[0,2].plot(ice_r * 1e6, z, label="ice", color=ice_c)
        ax[0,2].plot(liq_r * 1e6, z, label="liquid", color=liq_c)
        ax[0,3].plot(T, z, color=ice_c)
        ax[0,4].plot([rh_to_rh_i(RH_val, T_val) for (RH_val, T_val) in zip(RH, T) ], z, label='ice', color=ice_c)
        ax[0,4].plot(RH, z, label='liq', color=liq_c)

    ax[0,0].set_xlabel('LWC [g/m^3]')
    ax[0,0].set_ylabel('z [m]')
    ax[0,0].legend()
    ax[0,1].set_xlabel('concentration [1/mg]')
    ax[0,2].set_xlabel("average radius [um]")
    ax[0,3].set_xlabel("temperature [K]")
    ax[0,4].set_xlabel("RH")

    fig.tight_layout(rect=(0, 0.10, 1, 0.97))
    out_png = "test_sstp_WBF_mixing.png"
    plt.suptitle('WBF process with {sstp_cond} substeps')
    plt.savefig(out_png, dpi=200)

    return fig

make_figure('{"polluted": {"kappa": 0.61, "rd_insol" : 0.0, "mean_r": [0.029e-6, 0.071e-6], "gstdev": [1.36, 1.57], "n_tot": [160.0e6, 380.0e6]},' \
            '"INP": {"kappa": 0.61, "rd_insol" : 0.5e-6, "mean_r": [0.029e-6], "gstdev": [1.36], "n_tot": [10.0e6]}}')
plt.show()
