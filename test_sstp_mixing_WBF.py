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
from libcloudphxx import common
plt.rcParams.update({'font.size': 14})

sstp = 90
w_list = [1.]#[1., 2.5, 5.]
z_max_list = [3000]
outfile = f"test_WBF.nc"

polluted = '{"polluted": {"kappa": 0.61, "rd_insol" : 0.0, "mean_r": [0.029e-6, 0.071e-6], "gstdev": [1.36, 1.57], "n_tot": [160.0e6, 380.0e6]},' \
                '"INP": {"kappa": 0.61, "rd_insol" : 0.5e-6, "mean_r": [0.029e-6], "gstdev": [1.36], "n_tot": [10.0e6]}}' # low concentration of INPs

#pristine = '{"pristine": {"kappa": 0.61, "rd_insol": 0.0, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125.0e6, 65.0e6]},' \
#                '"INP": {"kappa": 0.61, "rd_insol" : 0.5e-6, "mean_r": [0.029e-6], "gstdev": [1.36], "n_tot": [10.0e6]}}' # low concentration of INPs

pristine = '{"pristine": {"kappa": 0.61, "rd_insol": 0.0, "mean_r": [0.011e-6], "gstdev": [1.2], "n_tot": [125.0e6]}}'


aerosol_list = [pristine]#, polluted]

monodisperse = {
    "ammonium_sulfate": {
        "kappa": 0.61, 
        "rd_insol": 0.0, 
        "bins": {1e-6: [30.0, 15]}
    }
}

def run_scheme(mixing, adaptive, aerosol, w_max, z_max):
    args = dict(
        p_0=100000,
        RH_0=0.7,
        T_0=283,
        aerosol = aerosol,
        #dry_sizes = monodisperse,
        sd_conc=100,
        dt=9,
        z_max=None,
        w = lambda t: w_max if t <= z_max/w_max else -w_max,
        t = z_max / w_max,
        outfile=outfile,
        outfreq=10,
        scheme="lgrngn",
        out_bin='{"liq": {"rght": 1, "moms": [0,1,2,3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0.5e-20},' \
                '"ice": {"rght": 1, "moms": [0,1,2,3], "drwt": "ice_a", "nbin": 1, "lnli": "lin", "left": 0.5e-20}}',
        sstp_cond = sstp,
        adaptive_sstp_cond = adaptive,
        sstp_cond_mix   = mixing,
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
        z = np.array(f.variables['z'][:]).squeeze()
        rv = np.array(f.variables['r_v'][:]).squeeze()
        RH = np.array(f.variables['RH'][:]).squeeze()
        T = np.array(f.variables['T'][:]).squeeze()
        act_m0 = np.array(f.variables['act_m0'][:]).squeeze()
        act_m1 = np.array(f.variables['act_m1'][:]).squeeze()
        act_m2 = np.array(f.variables['act_m2'][:]).squeeze()
        ice_m0 = np.array(f.variables['ice_m0'][:]).squeeze()
        ice_m1 = np.array(f.variables['ice_m1'][:]).squeeze()
        ice_m2 = np.array(f.variables['ice_m2'][:]).squeeze()
        ice_mix_ratio = np.array(f.variables['ice_mix_ratio'][:]).squeeze()
        liq_mix_ratio = np.array(f.variables['liq_m3'][:]).squeeze() * 4/3 * np.pi * common.rho_w
        ice_conc = np.array(f.variables['ice_m0'][:]).squeeze()  # 1/kg
        act_conc = np.array(f.variables['act_m0'][:]).squeeze()
        ice_r = np.where(ice_conc > 0, np.array(f.variables['ice_m1'][:]).squeeze() / np.array(f.variables['ice_m0'][:]).squeeze(), 0)
        act_r = np.where(act_conc > 0, np.array(f.variables['act_m1'][:]).squeeze() / np.array(f.variables['act_m0'][:]).squeeze(), 0)
        std_dev_liq = np.sqrt(np.where(act_m0 > 0, 
                           act_m2 / act_m0 - (act_m1 / act_m0)**2, 
                           0))
        std_dev_ice = np.sqrt(np.where(ice_m0 > 0, 
                           ice_m2 / ice_m0 - (ice_m1 / ice_m0)**2, 
                           0))
        sd_conc_liq = np.array(f.variables['sd_conc_liq'][:]).squeeze()
        sd_conc_ice = np.array(f.variables['sd_conc_ice'][:]).squeeze()
    #os.remove(outfile)
    return z/1000, rv*1e3, RH, T, ice_mix_ratio*1e3, liq_mix_ratio*1e3, ice_conc/1e6, act_conc/1e6, ice_r*1e6, act_r*1e6, std_dev_liq*1e6, std_dev_ice*1e6, sd_conc_liq, sd_conc_ice



def make_figure(aerosol, w_max, z_max):

    fig, ax = plt.subplots(2, 4, figsize=(15, 10.0), sharey=True, squeeze=True)

    for (mixing, adaptive) in [(False, False), (True, False)]:#, (True, True)]:#[(True, False), (True, True), (False, True)]:

        s = ':'
        lw = 2
        if mixing:
            l = 'per-particle coupled'
            c = 'blue'
        else:
            l = 'per-particle uncoupled'
            c = 'violet'

        z, rv, RH, T, ice_mix_ratio, liq_mix_ratio, ice_conc, act_conc, ice_r, act_r, std_dev_liq, std_dev_ice, sd_conc_liq, sd_conc_ice = run_scheme(mixing, adaptive, aerosol, w_max, z_max)            
        r_tot = rv + liq_mix_ratio + ice_mix_ratio

        #print('SD ice: ', sd_conc_ice)
        #print('SD  liq: ', sd_conc_liq)

        ax[1,0].plot(rv, z, color=c, label=l, linestyle=s, linewidth = lw)
        ax[0,0].plot(liq_mix_ratio, z, color=c, label=l, linestyle=s, linewidth = lw)
        ax[1,2].plot(T, z, color=c, label=l, linestyle=s, linewidth = lw)
        ax[0,1].plot(act_conc, z, color=c, label=l, linestyle=s, linewidth = lw)
        ax[1,3].plot(ice_mix_ratio, z, color=c, linestyle=s, linewidth = lw)
        ax[0,2].plot(act_r, z, color=c, linestyle=s, linewidth = lw)
        ax[0,3].plot(RH, z, color=c, label=l, linestyle=s, linewidth = lw)
        ax[1,1].plot(np.round(r_tot,1), z, color=c, label=l, linestyle=s, linewidth = lw)
        
    ax[0,0].set_ylabel('z [km]')
    ax[1,0].set_ylabel('z [km]')
    ax[0,0].set_xlabel('liquid mixing ratio [g/kg]')
    ax[1,0].set_xlabel('rv [g/kg]')
    ax[0,1].set_xlabel('liquid conc. [1/mg]')
    ax[1,2].set_xlabel('T [K]')
    ax[0,2].set_xlabel(f'liquid mean radius [$\mu$m]')
    ax[1,3].set_xlabel(f'ice mix ratio [g/kg]')
    ax[0,3].set_xlabel(f'RH')
    ax[1,1].set_xlabel(f'r_tot [g/kg]')

    handles, labels = ax[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center right", bbox_to_anchor=(0.21, 0.5))
    fig.tight_layout(rect=[0.2, 0, 1, 0.95])
    aerosol_str = "pristine" if aerosol==pristine else "polluted"
    out_png = "plots/outputs/coupled_uncoupled/test_mixing_w_"+str(w_max)+".svg"
    plt.suptitle('w = '+str(w_max)+' m/s, '+aerosol_str)
    plt.savefig(out_png, dpi=200)

    return fig

for w_max,z_max in zip(w_list, z_max_list):
    for aerosol in aerosol_list:
        make_figure(aerosol, w_max, z_max)
plt.show()