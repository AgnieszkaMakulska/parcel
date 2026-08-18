import coated_dust as cd
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8')
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
plt.rcParams.update({
    'font.size': 16,
    'axes.labelsize': 16,
    'axes.titlesize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16
})

aerosol_str = "polluted"

rd = 3e-6
epsilon = 0.1
z_distr = 400

mix = cd.mixed_aerosol(aerosol_str, epsilon, rd)
sol = cd.soluble_aerosol(aerosol_str)

out_png = "plots/rd_insol/spectrum_" + aerosol_str + ".pdf"
fig, ax = plt.subplots(1, 2, figsize=(12.0, 5.0), sharey=True, squeeze=False)
ax = ax.flatten()

for aerosol in [sol, mix]:
    outfile = "sol.nc" if aerosol == sol else "mix.nc"
    cd.run_scheme(aerosol, outfile, outfreq = z_distr, spec = True)
    distr1, radii, bin_widths, initial_distr, init_radii, init_bin_widths = cd.read_distr(outfile, z_distr)
    os.remove(outfile)

    l = "ammonium sulfate" if aerosol==sol else "ammonium sulfate + dust"
    c = "tab:blue" if aerosol==sol else "tab:orange"
    ax[0].bar(init_radii, initial_distr, color=c, edgecolor=c, width=init_bin_widths, alpha=0.4, label = l, linewidth=2)
    ax[1].bar(radii, distr1, color=c, edgecolor=c, width=bin_widths, alpha=0.4, label = l, linewidth=2)

ax[0].set_title('Initial size distribution')
ax[1].set_title('Size distribution at '+str(z_distr)+' m')
ax[0].set_xscale('log')
ax[1].set_xscale('log')
ax[0].set_yscale('log')
ax[1].set_yscale('log')
ax[0].set_ylabel('droplet concentration [1/mg]')        
ax[0].set_xlabel('droplet radius [$\\mu$m]')
ax[1].set_xlabel('droplet radius [$\\mu$m]')

if aerosol_str == "pristine":
    ax[1].set_xlim(12,19)
elif aerosol_str == "polluted":
    ax[1].set_xlim(6,19)

handles, labels = ax[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0.89),
            ncol=2, frameon=False)
plt.tight_layout(rect=[0, 0, 1, 0.90])

xmin, xmax = ax[1].get_xlim()
ticks = np.arange(np.ceil(xmin), np.floor(xmax) + 1)
ax[1].xaxis.set_major_locator(FixedLocator(ticks))
ax[1].xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
ax[1].xaxis.set_minor_locator(FixedLocator([]))
ax[1].xaxis.set_minor_formatter(NullFormatter())
ax[1].xaxis.get_offset_text().set_visible(False)
ax[1].tick_params(axis='x')

plt.savefig(out_png, dpi=200)