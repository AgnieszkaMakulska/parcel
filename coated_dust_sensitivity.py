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

from coated_dust import aerosol_spec, run_scheme, read_profiles

aerosol_str = "pristine"
out_png = "plots/rd_insol/sensitivity_" + aerosol_str + ".pdf"


rd_sol_list = np.linspace(0.0, 0.5, 10) * 1e-6
rd_insol_list = np.linspace(0.0, 10, 10) * 1e-6

lwc_list = []
nc_list = []
rc_list = []
a_stdev_list = []
x_coords = []
y_coords = []

for rd_sol in rd_sol_list:
    for rd_insol in rd_insol_list:

        if rd_insol == 0.0:
            aerosol = aerosol_spec(aerosol_str, 1.0)
        else:
            rd = np.cbrt(rd_insol**3 + rd_sol**3)
            epsilon = rd_sol**3 / rd**3
            aerosol = aerosol_spec(aerosol_str, epsilon, rd)
    
        outfile = outfile = str(rd_insol)+str(rd_sol)+".nc"
        run_scheme(aerosol, outfile, outfreq = 400)
        z, liq_mix_ratio, conc, mean_r, std_dev_area = read_profiles(outfile)
        os.remove(outfile)

        lwc_list.append(liq_mix_ratio[-1])
        nc_list.append(conc[-1])
        rc_list.append(mean_r[-1])
        a_stdev_list.append(std_dev_area[-1])
        x_coords.append(rd_sol)
        y_coords.append(rd_insol)

n_sol = len(rd_sol_list)
n_insol = len(rd_insol_list)

def to_grid(flat_list):
    arr = np.array(flat_list).reshape(n_sol, n_insol).T
    return arr


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

    cbar = fig.colorbar(im, ax=ax[i])
    ax[i].set_title(titles[i])

for i in (2,3):
    ax[i].set_xlabel("$r_{sol}$ [$\\mu$m]")
for i in (0,2):
    ax[i].set_ylabel("$r_{insol}$ [$\\mu$m]")

plt.tight_layout()
plt.savefig(out_png, dpi=200, bbox_inches='tight')

    



