import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")

import numpy as np
from parcel import parcel
from scipy.io import netcdf
import matplotlib.pyplot as plt
import os
from libcloudphxx import common
plt.style.use('seaborn-v0_8')
plt.rcParams.update({
    'font.size': 16,
    'axes.labelsize': 16,
    'axes.titlesize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16
})

w = 1.
z_max = 400

# polluted = '{"polluted": {"kappa": 1.28, "sol_frac": 1.0, "mean_r": [0.029e-6, 0.071e-6], "gstdev": [1.36, 1.57], "n_tot": [160.0e6, 380.0e6]}}'
# pristine = '{"pristine": {"kappa": 1.28, "sol_frac": 1.0, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125.0e6, 65.0e6]}}'


def run_scheme(epsilon, rd, outfile):
    
    if epsilon == 1.0:
        aerosol = '{"pristine": {"kappa": 1.28, "sol_frac": 1.0, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125.0e6, 65.0e6]}}'
    else:
        aerosol = f'{{"soluble":{{"kappa": 1.28, "sol_frac": 1.0, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125.0e6, 65.0e6]}}, \
                "mixed": {{"kappa": 1.28, "sol_frac": {epsilon}, "mean_r": [{rd}], "gstdev": [1.4], "n_tot": [1.0e6]}} }}'

    args = dict(
        p_0=90000,
        RH_0=0.97,
        T_0=283,
        aerosol = aerosol,
        w = w,
        sd_conc=1000,
        #sd_const_multi=1000000,
        #n_sd_max=1e7,
        dt=1,
        z_max = z_max,
        outfile=outfile,
        outfreq=1,
        scheme="lgrngn",
        out_bin='{"liq": {"rght": 1, "moms": [0,1,2,3,4], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0.5e-20}}',
        sstp_cond = 10,
        ice_switch = False,
        ice_nucl = False,
        depo = False,
        backend = "gpu"
    )
    parcel(**args)


def read(outfile):
    with netcdf.netcdf_file(outfile, 'r') as f:
        z = np.array(f.variables['z'][:]).squeeze()
        act_m0 = np.array(f.variables['act_m0'][:]).squeeze()
        act_m1 = np.array(f.variables['act_m1'][:]).squeeze()
        act_m2 = np.array(f.variables['act_m2'][:]).squeeze()
        act_m4 = np.array(f.variables['act_m4'][:]).squeeze()
        liq_m0 = np.array(f.variables['liq_m0'][:]).squeeze()
        liq_m1 = np.array(f.variables['liq_m1'][:]).squeeze()
        liq_m2 = np.array(f.variables['liq_m2'][:]).squeeze()
        liq_m3 = np.array(f.variables['liq_m3'][:]).squeeze()
        liq_m4 = np.array(f.variables['liq_m4'][:]).squeeze()
        liq_mix_ratio = liq_m3 * 4/3 * np.pi * common.rho_w
        conc = act_m0
        mean_r = np.where(act_m0 > 0, act_m1 / act_m0, 0)
        std_dev_area = np.sqrt(np.where(act_m0 > 0, 
                        act_m4 / act_m0 - (act_m2 / act_m0)**2, 
                        0)) * 4 * np.pi
    return z, liq_mix_ratio*1e3, conc/1e6, mean_r*1e6, std_dev_area*1e12



def make_profiles():

    rd = 1e-6
    epsilon_list = [0.0, 1e-3, 1e-2, 1e-1, 1.]

    out_png = "plots/rd_insol/different_epsilon.pdf"

    fig, ax = plt.subplots(2, 2, figsize=(8.0, 9.0), sharey=True, squeeze=False)
    ax = ax.flatten()
    
    for epsilon in epsilon_list:

        l = '$\\epsilon$ = ' + str(round(epsilon, 3))

        outfile = str(epsilon)+".nc"
        run_scheme(epsilon, rd, outfile)
        z, liq_mix_ratio, conc, mean_r, std_dev_area = read(outfile)
        os.remove(outfile)
     
        ax[0].plot(liq_mix_ratio, z, label=l)
        ax[1].plot(conc, z, label=l)
        ax[2].plot(mean_r, z, label=l)
        ax[3].plot(std_dev_area, z, label=l)

    ax[0].set_ylabel('z [m]')
    ax[2].set_ylabel('z [m]')
    ax[0].set_xlabel('liquid mix. ratio [g/kg]')
    ax[1].set_xlabel('droplet concentration [1/mg]')
    ax[2].set_xlabel(f'droplet mean radius [$\mu$m]')
    ax[3].set_xlabel(f'std. dev. of droplet area [$\mu$m$^2$]')
    ax[0].set_xlim(0.1,0.6)
    ax[1].set_xlim(40,65)
    ax[2].set_xlim(5,15)
    ax[3].set_xlim(40,90)

    handles, labels = ax[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0.89),
               ncol=2, frameon=False)
    plt.tight_layout(rect=[0, 0, 1, 0.9])

    plt.savefig(out_png, dpi=200, bbox_inches='tight')

make_profiles()