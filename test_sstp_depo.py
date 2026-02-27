"""
Checking if depositional growth is sensitive to the number of substeps
"""

import sys, os
sys.path.insert(0, "../")
sys.path.insert(0, "./")

import numpy as np
from parcel import parcel
from scipy.io import netcdf
import matplotlib.pyplot as plt
from functions import rh_to_rh_i

sstp_list = [1, 5, ]
#sstp_list = [1]

aerosol = '{"polluted": {"kappa": 0.61, "rd_insol" : 0.5e-6, "mean_r": [0.029e-6, 0.071e-6], "gstdev": [1.36, 1.57], "n_tot": [160.0e6, 380.0e6]}}'

def run_scheme(outfile, sstp):
    args = dict(
        p_0=100000,
        RH_0=0.9,
        T_0=260,
        aerosol = aerosol,
        sd_conc=100,
        dt=1,
        z_max=1000,
        w = 1.,
        outfile=outfile,
        outfreq=10,
        scheme="lgrngn",
        out_bin='{"liq": {"rght": 1, "moms": [0,1,3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0.5e-10},' \
                '"ice": {"rght": 1, "moms": [0,1,3], "drwt": "ice_a", "nbin": 1, "lnli": "lin", "left": 0.5e-10}}',
        sstp_cond = sstp,
        adaptive_sstp_cond = False,
        sstp_cond_mix   = False,
        exact_sstp_cond = True,
        aerosol_independent_of_rhod=True, 
        backend="OpenMP",
        ice_switch = True,
        ice_nucl = True,
        time_dep_ice_nucl = True,
        depo = True
    )

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



def make_figure():

    fig, ax = plt.subplots(2, 4, figsize=(12.0, 12.0), sharey=True, squeeze=False)


    for sstp in sstp_list:

        outfile = f"test_dep_sstp_{sstp}.nc"
        RH, T, z, ice_mix_ratio, liq_mix_ratio, liq_conc, ice_conc, liq_r, ice_r = run_scheme(outfile, sstp)

        ax[0,0].plot(ice_mix_ratio * 1e3, z, label = str(sstp))
        ax[0,1].plot(liq_mix_ratio * 1e3, z)
        ax[0,2].plot(ice_conc / 1e6, z)
        ax[0,3].plot(liq_conc / 1e6, z)
        ax[1,0].plot(ice_r * 1e6, z)
        ax[1,1].plot(liq_r * 1e6, z)
        ax[1,2].plot(T, z)
        ax[1,3].plot(RH, z)

    ax[0,0].legend(loc='center')
    ax[0,0].set_xlabel('IWC [g/m^3]')
    ax[0,1].set_xlabel('LWC [g/m^3]')
    ax[0,2].set_xlabel('ice concentration [1/mg]')
    ax[0,3].set_xlabel('liq concentration [1/mg]')
    ax[1,0].set_xlabel("avg ice radius [um]")
    ax[1,1].set_xlabel("avg liq radius [um]")
    ax[1,2].set_xlabel("temperature [K]")
    ax[1,3].set_xlabel("RH")

    fig.tight_layout(rect=(0, 0.10, 1, 0.97))
    out_png = "test_depo_sstp.png"
    plt.savefig(out_png, dpi=200)

    return fig


make_figure()

plt.show()
