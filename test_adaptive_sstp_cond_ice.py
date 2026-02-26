"""
Run adaptive substepping with ice
"""

import sys, os
sys.path.insert(0, "../")
sys.path.insert(0, "./")

import numpy as np
from parcel import parcel
from scipy.io import netcdf
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import matplotlib.colors as mcolors
from pathlib import Path
from typing import List

sstp_cond_max = 10
z_max = 3000

def run_scheme(w_max, outfile, *, sstp_cond=sstp_cond_max):
    args = dict(
        p_0=100000,
        RH_0=0.9,
        T_0=273,
        aerosol = None,
        sd_conc=100,
        dt=1,
        z_max=None,
        w = w_max,
        outfile=outfile,
        outfreq=10,
        scheme="lgrngn",
        out_bin='{"liq": {"rght": 1, "moms": [0,1,3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0.5e-6},' \
                '"ice": {"rght": 1, "moms": [0,1,3], "drwt": "ice_a", "nbin": 1, "lnli": "lin", "left": 0.5e-6}}',
        sstp_cond=sstp_cond,
        adaptive_sstp_cond=False,
        sstp_cond_mix   = True,
        exact_sstp_cond = False,
        aerosol_independent_of_rhod=True, 
        backend="OpenMP",
        ice_switch = True,
        ice_nucl = True,
        time_dep_ice_nucl = False,
        depo = True,
        wait = 0
    )

    if hasattr(run_scheme, "aerosol"):
        args["aerosol"] = run_scheme.aerosol
    if hasattr(run_scheme, "sstp_cond_adapt_drw2_eps"):
        args["sstp_cond_adapt_drw2_eps"] = float(run_scheme.sstp_cond_adapt_drw2_eps)
    if hasattr(run_scheme, "sstp_cond_adapt_drw2_max"):
        args["sstp_cond_adapt_drw2_max"] = float(run_scheme.sstp_cond_adapt_drw2_max)
    if hasattr(run_scheme, "sstp_cond_act"):
        args["sstp_cond_act"] = int(run_scheme.sstp_cond_act)      


    #args["t"] = 2. * z_max / w_max
    args["t"] = z_max / w_max
    args["outfreq"] = 1
    print("t: ", args["t"])
    parcel(**args)

    with netcdf.netcdf_file(outfile, 'r') as f:
        z = np.array(f.variables['z'][:])
        RH = np.array(f.variables['RH'][:])
        T = np.array(f.variables['T'][:])
        sstp_cond_mean = np.array(f.variables['sstp_cond_mean'][:]) if 'sstp_cond_mean' in f.variables else None
        if sstp_cond_mean is not None:
            sstp_cond_mean[0] = sstp_cond_mean[1] # at t=0 sstp_cond_mean=0, because its set only during the firs step (?)
        ice_mix_ratio = np.array(f.variables['ice_mix_ratio'][:])
        liq_mix_ratio = np.array(f.variables['liq_m3'][:]) *4/3 * np.pi * 997 #multiply by density of water
        liq_conc = np.array(f.variables['liq_m0'][:])  # 1/kg
        ice_conc = np.array(f.variables['ice_m0'][:])  # 1/kg
        liq_r = np.where(liq_conc > 0, np.array(f.variables['liq_m1'][:]) / np.array(f.variables['liq_m0'][:]), 0)
        ice_r = np.where(ice_conc > 0, np.array(f.variables['ice_m1'][:]) / np.array(f.variables['ice_m0'][:]), 0)
    return RH, T, z, sstp_cond_mean, ice_mix_ratio, liq_mix_ratio, liq_conc, ice_conc, liq_r, ice_r


# baseline - basically no adaptation, very relaxed conditions
baseline = dict(
    eps=1e6, #1e-1,
    max=1e6, #100,
    act=1,  # 1 means disabled
)


def make_figure(aerosol_name, aerosol):
    run_scheme.aerosol = aerosol
    w_max = 1.0
    eps = 1e-2
    fig, ax = plt.subplots(1, 5, figsize=(15.0, 15.0), sharey=True, squeeze=False)

    generated_nc_files: List[str] = []

    cmap_dt = "gnuplot"
    norm_dt = mcolors.Normalize(vmin=1, vmax=sstp_cond_max)


    run_scheme.sstp_cond_adapt_drw2_eps = eps
    run_scheme.sstp_cond_adapt_drw2_max = baseline["max"]
    run_scheme.sstp_cond_act = baseline["act"]

    outfile = f"test_adaptive_sstp_cond_{aerosol_name}_w{w_max:g}_eps{eps:.0e}_adapt1.nc"
    RH, T, z, sstp_cond_mean, ice_mix_ratio, liq_mix_ratio, liq_conc, ice_conc, liq_r, ice_r = run_scheme(w_max, outfile)
    generated_nc_files.append(outfile)


    ax[0,0].plot(ice_mix_ratio, z, label="ice")
    ax[0,0].plot(liq_mix_ratio, z, label="liq")
    ax[0,0].legend()
    ax[0,0].set_title(f"eps={eps:.0e}")
    ax[0,0].set_ylabel(f"w_max={w_max:g}\nHeight [m]")
    ax[0,0].set_xlabel('LWC')

    ax[0,1].plot(ice_conc / 1e6, z, label="ice")
    ax[0,1].plot(liq_conc / 1e6, z, label="liquid")
    ax[0,1].legend()
    ax[0,1].set_xlabel("number conc. [1/mg]")

    ax[0,2].plot(ice_r * 1e6, z, label="ice")
    ax[0,2].plot(liq_r * 1e6, z, label="liquid")
    ax[0,2].legend()
    ax[0,2].set_xlabel("average radius [um]")

    ax[0,3].plot(T, z)
    ax[0,3].set_xlabel("Temperature [K]")

    from functions import rh_to_rh_i
    ax[0,4].plot([rh_to_rh_i(RH_val, T_val) for (RH_val, T_val) in zip(RH, T) ], z, label='ice')
    ax[0,4].plot(RH, z, label='liq')
    ax[0,4].set_xlabel("RH")
    ax[0,4].legend()


    fig.suptitle("Adaptive substepping, " + aerosol_name)
    fig.tight_layout(rect=(0, 0.10, 1, 0.97))

    out_png = "test_adaptive_sstp_cond_"+aerosol_name+".png"
    plt.savefig(out_png, dpi=200)

    return fig

#make_figure('pristine', '{"DYCOMS": {"kappa": 0.61, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125.0e6, 65.0e6]}}', 100)
make_figure('polluted', '{"polluted": {"kappa": 0.61, "rd_insol" : 0.0, "mean_r": [0.029e-6, 0.071e-6], "gstdev": [1.36, 1.57], "n_tot": [160.0e6, 380.0e6]},' \
            '"INP": {"kappa": 0.61, "rd_insol" : 0.5e-6, "mean_r": [0.029e-6], "gstdev": [1.36], "n_tot": [160.0e6]}}')
plt.show()
