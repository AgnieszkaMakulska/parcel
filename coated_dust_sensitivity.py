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
    return liq_mix_ratio[-1]*1e3, conc[-1]/1e6, mean_r[-1]*1e6, std_dev_area[-1]*1e12


rd_sol_list = np.linspace(0.01, 0.1, 5) * 1e-6
rd_insol_list = np.linspace(0.1, 10, 5) * 1e-6

lwc_list = []
nc_list = []
rc_list = []
a_stdev_list = []
x_coords = []
y_coords = []

for rd_sol in rd_sol_list:
    for rd_insol in rd_insol_list:
        if rd_insol == 0.0:
            rd = rd_sol
            epsilon = 1.0
        else:
            rd = np.cbrt(rd_insol**3 + rd_sol**3)
            epsilon = rd_sol**3 / rd**3
    
        outfile = "test.nc"
        run_scheme(epsilon, rd, outfile)
        liq_mix_ratio, conc, mean_r, std_dev_area = read(outfile)
        os.remove(outfile)

        lwc_list.append(liq_mix_ratio)
        nc_list.append(conc)
        rc_list.append(mean_r)
        a_stdev_list.append(std_dev_area)
        x_coords.append(rd_sol)
        y_coords.append(rd_insol)

n_sol = len(rd_sol_list)
n_insol = len(rd_insol_list)

def to_grid(flat_list):
    arr = np.array(flat_list).reshape(n_sol, n_insol).T  # -> shape (n_insol, n_sol)
    return arr

# print(to_grid(x_coords))
# print(to_grid(y_coords))

datasets = [to_grid(lwc_list), to_grid(nc_list), to_grid(rc_list), to_grid(a_stdev_list)]
titles = [
    'liquid mix. ratio [g/kg]',
    'droplet concentration [1/mg]',
    f'droplet mean radius [$\mu$m]',
    f'std. dev. of droplet area [$\mu$m$^2$]'
]

fig, ax = plt.subplots(2, 2, figsize=(14.0, 12.0), sharey=False, squeeze=True)
ax = ax.flatten()

for i in range(4):
    X, Y = np.meshgrid(range(n_sol+1), range(n_insol+1))
    im = ax[i].pcolormesh(X, Y, datasets[i], cmap='viridis')

    ax[i].set_xticks(range(n_sol))
    ax[i].set_xticklabels([f"{v*1e6:.2f}" for v in rd_sol_list], rotation=45)
    ax[i].set_yticks(range(n_insol))
    ax[i].set_yticklabels([f"{v*1e6:.1f}" for v in rd_insol_list])

    ax[i].set_xlabel("$r_{sol}$ [$\\mu$m]")
    ax[i].set_ylabel("$r_{insol}$ [$\\mu$m]")
    cbar = fig.colorbar(im, ax=ax[i])
    cbar.set_label(titles[i])

plt.tight_layout()
out_png = "plots/rd_insol/sensitivity.pdf"
plt.savefig(out_png, dpi=200, bbox_inches='tight')

    



