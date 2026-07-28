"""
Checking how insoluble component impacts condensation
"""

import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")

import numpy as np
from parcel import parcel
from scipy.io import netcdf
import matplotlib.pyplot as plt
from libcloudphxx import common
plt.rcParams.update({'font.size': 16})
from matplotlib.ticker import LogLocator, FuncFormatter, NullFormatter
import json

w = 1.
z_max = 500

# polluted = '{"polluted": {"kappa": 1.28, "sol_frac": 1.0, "mean_r": [0.029e-6, 0.071e-6], "gstdev": [1.36, 1.57], "n_tot": [160.0e6, 380.0e6]}}'
# pristine = '{"pristine": {"kappa": 1.28, "sol_frac": 1.0, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125.0e6, 65.0e6]}}'


def run_scheme(epsilon, rd, outfile):

    mix = f'{{"soluble":{{"kappa": 1.28, "sol_frac": 1.0, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125.0e6, 65.0e6]}}, \
                "mixed": {{"kappa": 1.28, "sol_frac": {epsilon}, "mean_r": [{rd}], "gstdev": [1.4], "n_tot": [1.0e6]}} }}'

    args = dict(
        p_0=90000,
        RH_0=0.97,
        T_0=283,
        aerosol = mix,
        w = w,
        sd_conc=100,
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
        time_dep_ice_nucl = True,
        depo = False,
        backend = "OpenMP"
    )
    parcel(**args)


def read(outfile):
    with netcdf.netcdf_file(outfile, 'r') as f:
        act_m0 = np.array(f.variables['act_m0'][-1]).squeeze()
        act_m1 = np.array(f.variables['act_m1'][-1]).squeeze()
        act_m2 = np.array(f.variables['act_m2'][-1]).squeeze()
        act_m4 = np.array(f.variables['act_m4'][-1]).squeeze()
        liq_m0 = np.array(f.variables['liq_m0'][-1]).squeeze()
        liq_m1 = np.array(f.variables['liq_m1'][-1]).squeeze()
        liq_m2 = np.array(f.variables['liq_m2'][-1]).squeeze()
        liq_m3 = np.array(f.variables['liq_m3'][-1]).squeeze()
        liq_m4 = np.array(f.variables['liq_m4'][-1]).squeeze()
        liq_mix_ratio = liq_m3 * 4/3 * np.pi * common.rho_w
        conc = act_m0
        mean_r = np.where(act_m0 > 0, act_m1 / act_m0, 0)
        std_dev_area = np.sqrt(np.where(act_m0 > 0, 
                        act_m4 / act_m0 - (act_m2 / act_m0)**2, 
                        0)) * 4 * np.pi
    return liq_mix_ratio*1e3, conc/1e6, mean_r*1e6, std_dev_area*1e12



rd_sol_list = np.array([0.01, 0.04, 0.08]) * 1e-6
rd_insol_list = np.array([0.3, 1., 3., 5.]) * 1e-6

lwc_list = []
nc_list = []
rc_list = []
a_stdev_list = []
x_coords = []
y_coords = []

for rd_sol in rd_sol_list:
    for rd_insol in rd_insol_list:
        rd = np.cbrt(rd_insol**3 + rd_sol**3)
        epsilon = rd_sol**3 / rd**3
        outfile = "test.nc"

        run_scheme(epsilon, rd, outfile)
        liq_mix_ratio, conc, mean_r, std_dev_area = read(outfile)

        lwc_list.append(liq_mix_ratio)
        nc_list.append(conc)
        rc_list.append(mean_r)
        a_stdev_list.append(std_dev_area)
        x_coords.append(rd_sol)
        y_coords.append(rd_insol)

datasets = [lwc_list, nc_list, rc_list, a_stdev_list]
titles = [
    "liq mix ratio [g/kg]",
    "act conc [1/cm^3]",
    "mean radius [um]",
    "area stdev [um^2]"
]

fig, ax = plt.subplots(2, 2, figsize=(16.0, 15.0), sharey=False, squeeze=True)
ax = ax.flatten()

for i in range(4):
    sc = ax[i].scatter(
        np.array(x_coords) * 1e6,
        np.array(y_coords) * 1e6,
        c=datasets[i],
        cmap='viridis',
        alpha=0.8
    )

    ax[i].set_xlabel("rd_sol")
    ax[i].set_ylabel("rd_insol")
    cbar = fig.colorbar(sc, ax=ax[i])
    cbar.set_label(titles[i])

plt.tight_layout()
plt.show()


    



