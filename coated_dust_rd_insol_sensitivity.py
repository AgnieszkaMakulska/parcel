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

aerosol_str = "polluted"
out_png = "plots/rd_insol/different_rd_insol_" + aerosol_str +".pdf"

fig, ax = plt.subplots(2, 2, figsize=(8.0, 9.0), sharey=True, squeeze=False)
ax = ax.flatten()

rd_insol_list = np.array([0.0, 0.1, 1., 8.]) * 1e-6
for rd_insol in rd_insol_list:

    if rd_insol == 0.0:
        l = 'no dust'
        aerosol = aerosol_spec(aerosol_str, 1.0)
    
    else:
        epsilon = 0.01
        l = '$r_{insol}$ = ' + str(round(rd_insol*1e6, 1)) + ' $\\mu$m'
        rd = np.cbrt( rd_insol**3 / (1-epsilon) )
        aerosol = aerosol_spec(aerosol_str, epsilon, rd)
        
    outfile = str(rd_insol)+".nc"
    run_scheme(aerosol, outfile, outfreq = 5)
    z, liq_mix_ratio, conc, mean_r, std_dev_area = read_profiles(outfile)
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
if aerosol_str == "pristine":
    ax[1].set_xlim(40,65)
    ax[2].set_xlim(5,15)
    ax[3].set_xlim(40,90)
elif aerosol_str == "polluted":
    ax[1].set_xlim(250,410)
    ax[2].set_xlim(4,8)
    ax[3].set_xlim(25,38)


handles, labels = ax[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0.89),
            ncol=2, frameon=False)
plt.tight_layout(rect=[0, 0, 1, 0.9])

plt.savefig(out_png, dpi=200, bbox_inches='tight')
