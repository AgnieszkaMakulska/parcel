import coated_dust as cd
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8')
plt.rcParams.update({
    'font.size': 16,
    'axes.labelsize': 16,
    'axes.titlesize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16
})

aerosol_str = "pristine"

rd_list = np.logspace(-2, 0.7, 5) * 1e-6
eps_list = np.linspace(0.01, 0.9, 5)

lwc_list = []
nc_list = []
rc_list = []
r_stdev_list = []

for rd in rd_list:
    for epsilon in eps_list:
        print(rd,epsilon)

        aerosol = cd.mixed_aerosol(aerosol_str, epsilon, rd)
        outfile = str(rd)+str(epsilon)+".nc"
        cd.run_scheme(aerosol, outfile, outfreq = 400)
        z, rh, liq_mix_ratio, conc, mean_r, std_dev_r = cd.read_profiles(outfile)
        os.remove(outfile)

        lwc_list.append(liq_mix_ratio[-1])
        nc_list.append(conc[-1])
        rc_list.append(mean_r[-1])
        r_stdev_list.append(std_dev_r[-1])

n_rd = len(rd_list)
n_eps = len(eps_list)

def to_grid(flat_list):
    arr = np.array(flat_list).reshape(n_rd, n_eps).T
    return arr


datasets = [to_grid(lwc_list), to_grid(nc_list), to_grid(rc_list), to_grid(r_stdev_list)]
titles = [
    'liquid mix. ratio [g/kg]',
    'droplet concentration [1/mg]',
    'droplet mean radius [$\\mu$m]',
    'std. dev. of droplet radius [$\\mu$m]'
]

out_png = "plots/rd_insol/sensitivity_" + aerosol_str + ".pdf"
fig, ax = plt.subplots(2, 2, figsize=(14.0, 12.0), sharey=False, squeeze=True)
ax = ax.flatten()

for i in range(4):
    X, Y = np.meshgrid(range(n_rd+1), range(n_eps+1))
    im = ax[i].pcolormesh(X, Y, datasets[i], cmap='summer')

    ax[i].set_xticks(range(n_rd))
    ax[i].set_xticklabels([f"{v*1e6:.2f}" for v in rd_list], rotation=45)
    ax[i].set_yticks(range(n_eps))
    ax[i].set_yticklabels([f"{v:.2f}" for v in eps_list])

    cbar = fig.colorbar(im, ax=ax[i])
    ax[i].set_title(titles[i])

for i in range(4):
    ax[i].set_xlabel("$r_d$ [$\\mu$m]")
for i in range(4):
    ax[i].set_ylabel("$\\epsilon$")

plt.tight_layout()
plt.savefig(out_png, dpi=200, bbox_inches='tight')

    



