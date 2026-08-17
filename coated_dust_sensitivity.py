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


rd_list = np.linspace(0.01, 5, 10) * 1e-6
eps_list = np.linspace(0.0, 1., 10)

lwc_list = []
nc_list = []
rc_list = []
r_stdev_list = []
x_coords = []
y_coords = []

for rd in rd_list:
    for epsilon in eps_list:

        aerosol = aerosol_spec(aerosol_str, epsilon, rd)
    
        outfile = outfile = str(rd)+str(epsilon)+".nc"
        run_scheme(aerosol, outfile, outfreq = 400)
        z, liq_mix_ratio, conc, mean_r, std_dev_r = read_profiles(outfile)
        os.remove(outfile)

        lwc_list.append(liq_mix_ratio[-1])
        nc_list.append(conc[-1])
        rc_list.append(mean_r[-1])
        r_stdev_list.append(std_dev_r[-1])
        x_coords.append(rd)
        y_coords.append(epsilon)

n_rd = len(rd_list)
n_eps = len(eps_list)

def to_grid(flat_list):
    arr = np.array(flat_list).reshape(n_rd, n_eps).T
    return arr


datasets = [to_grid(lwc_list), to_grid(nc_list), to_grid(rc_list), to_grid(r_stdev_list)]
titles = [
    'liquid mix. ratio [g/kg]',
    'droplet concentration [1/mg]',
    f'droplet mean radius [$\\mu$m]',
    f'std. dev. of droplet radius [$\\mu$m]'
]

fig, ax = plt.subplots(2, 2, figsize=(14.0, 12.0), sharey=False, squeeze=True)
ax = ax.flatten()

for i in range(4):
    X, Y = np.meshgrid(range(n_rd+1), range(n_eps+1))
    im = ax[i].pcolormesh(X, Y, datasets[i], cmap='summer')

    ax[i].set_xticks(range(n_rd))
    ax[i].set_xticklabels([f"{v*1e6:.2f}" for v in rd_list], rotation=45)
    ax[i].set_yticks(range(n_eps))
    ax[i].set_yticklabels([f"{v:.1f}" for v in eps_list])

    cbar = fig.colorbar(im, ax=ax[i])
    ax[i].set_title(titles[i])

for i in (2,3):
    ax[i].set_xlabel("$r_d$ [$\\mu$m]")
for i in (0,2):
    ax[i].set_ylabel("$\epsilon$")

plt.tight_layout()
plt.savefig(out_png, dpi=200, bbox_inches='tight')

    



