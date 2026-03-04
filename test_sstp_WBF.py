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
from libcloudphxx import common
from matplotlib.ticker import MaxNLocator

sstp = 10
w_max = 1.
z_max = 1500.

aerosol = '{"polluted": {"kappa": 0.61, "rd_insol" : 0.0, "mean_r": [0.029e-6, 0.071e-6], "gstdev": [1.36, 1.57], "n_tot": [160.0e6, 380.0e6]},' \
                '"INP": {"kappa": 0.61, "rd_insol" : 0.5e-6, "mean_r": [0.029e-6], "gstdev": [1.36], "n_tot": [10.0e6]}}' # low concentration of INPs

def run_scheme(outfile, mixing, per_particle):
    args = dict(
        p_0=100000,
        RH_0=0.9,
        T_0=277,
        aerosol = aerosol,
        sd_conc=100,
        dt=1,
        # z_max=z_max,
        # w = 1.,
        z_max=None,
        w=lambda t: w_max * np.pi / 2. * np.sin(np.pi*t*w_max/z_max),
        t = z_max / w_max,
        outfile=outfile,
        outfreq=10,
        scheme="lgrngn",
        out_bin='{"liq": {"rght": 1, "moms": [0,1,3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0.5e-20},' \
                '"ice": {"rght": 1, "moms": [0,1,3], "drwt": "ice_a", "nbin": 1, "lnli": "lin", "left": 0.5e-20}}',
        sstp_cond = sstp,
        adaptive_sstp_cond = False,
        sstp_cond_mix   = mixing,
        exact_sstp_cond = per_particle,
        aerosol_independent_of_rhod=True, 
        backend="OpenMP",
        ice_switch = True,
        ice_nucl = True,
        time_dep_ice_nucl = False,
        depo = True
    )

    parcel(**args)

    with netcdf.netcdf_file(outfile, 'r') as f:
        z = np.array(f.variables['z'][:]).squeeze()
        RH = np.array(f.variables['RH'][:]).squeeze()
        T = np.array(f.variables['T'][:]).squeeze()
        rv = np.array(f.variables['r_v'][:]).squeeze()
        ice_mix_ratio = np.array(f.variables['ice_mix_ratio'][:]).squeeze()
        liq_mix_ratio = np.array(f.variables['liq_m3'][:]).squeeze() * 4/3 * np.pi * common.rho_w  #multiply by density of water
        ice_conc = np.array(f.variables['ice_m0'][:]).squeeze()  # 1/kg
        act_conc = np.array(f.variables['act_m0'][:]).squeeze()
        ice_r = np.where(ice_conc > 0, np.array(f.variables['ice_m1'][:]).squeeze() / np.array(f.variables['ice_m0'][:]).squeeze(), 0)
        act_r = np.where(act_conc > 0, np.array(f.variables['act_m1'][:]).squeeze() / np.array(f.variables['act_m0'][:]).squeeze(), 0)
    return RH, T, rv, z, ice_mix_ratio, liq_mix_ratio, ice_conc, act_conc, ice_r, act_r



def make_figure():

    fig, ax = plt.subplots(1, 5, figsize=(15.0, 8.0), sharey=True, squeeze=False)

    for (mixing, per_particle) in [(True, True), (False, True), (True, False)]:

        outfile = f"test_WBF.nc"
        RH, T, rv, z, ice_mix_ratio, liq_mix_ratio, ice_conc, act_conc, ice_r, act_r = run_scheme(outfile, mixing, per_particle)            

        if not per_particle:
            (ice_c, liq_c, ice_l, liq_l) = ("dimgrey", "dimgrey", "ice (per-cell)", "liquid (per-cell)")
            (c, l) = ("dimgrey", "per-cell")
            s = ':'
        if per_particle and mixing:
            (ice_c, liq_c, ice_l, liq_l) = ("skyblue", "coral", "ice (coupled)", "liquid (coupled)")
            (c, l) = ("violet", "coupled")
            s = '-'
        if per_particle and not mixing:
            (ice_c, liq_c, ice_l, liq_l) = ("steelblue", "sienna", "ice (uncoupled)", "liquid (uncoupled)")
            (c, l) = ("purple", "uncoupled")
            s = '-'
        
        ax[0,0].plot(ice_mix_ratio * 1e3, z, color=ice_c, label=ice_l, linestyle=s)
        ax[0,0].plot(liq_mix_ratio * 1e3, z, color=liq_c, label=liq_l, linestyle=s)
        ax[0,1].plot(ice_conc / 1e6, z, label=ice_l, color=ice_c, linestyle=s)
        ax[0,1].plot(act_conc / 1e6, z, label=liq_l, color=liq_c, linestyle=s)
        ax[0,2].plot(ice_r * 1e6, z, color=ice_c, linestyle=s)
        ax[0,2].plot(act_r * 1e6, z, color=liq_c, linestyle=s)
        ax[0,3].plot((rv + liq_mix_ratio + ice_mix_ratio) * 1e3, z, color=c, label=l, linestyle=s)
        ax[0,4].plot([rh_to_rh_i(RH_val, T_val) for (RH_val, T_val) in zip(RH, T) ], z, color=ice_c, linestyle=s)
        ax[0,4].plot(RH, z, color=liq_c, linestyle=s)

    ax[0,3].legend()
    handles, labels = ax[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.15, 0.2))

    ax[0,0].set_ylabel('z [m]')
    ax[0,0].set_xlabel('mix ratio [g/kg]')
    ax[0,1].set_xlabel('concentration [1/mg]')
    ax[0,2].set_xlabel("average radius [um]")
    ax[0,3].set_xlabel("total mix ratio [g/kg]")
    ax[0,4].set_xlabel("RH")

    ax[0,3].ticklabel_format(style='plain', useOffset=False)
    ax[0,3].xaxis.set_major_locator(MaxNLocator(4))

    fig.tight_layout(rect=(0, 0.10, 1, 0.97))
    out_png = "test_WBF_mixing.png"
    plt.suptitle('WBF process with ' + str(sstp) + ' substeps')
    plt.savefig(out_png, dpi=200)

    return fig

make_figure()
plt.show()