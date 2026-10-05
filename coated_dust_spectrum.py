import coated_dust as cd
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8-colorblind')
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter
plt.rcParams.update({
    'font.size': 16,
    'axes.labelsize': 16,
    'axes.titlesize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16
})

aerosol_str = "pristine"

epsilon = 0.05
zmax = 400

mix = cd.mixed_aerosol(aerosol_str, epsilon)
sol = cd.soluble_aerosol(aerosol_str)

sol_str = "sea salt" if aerosol_str == "pristine" else "ammonium sulfate"

out_png = "plots/rd_insol/spectrum_" + aerosol_str + ".pdf"
fig, ax = plt.subplots(2, 2, figsize=(12.0, 10.0), sharey=True, squeeze=False)
ax = ax.flatten()

for aerosol in [sol, mix]:
    outfile = "sol.nc" if aerosol == sol else "mix.nc"
    cd.run_scheme(aerosol, zmax, outfile, outfreq = zmax, spec = True)
    distr1, radii, bin_widths, initial_distr, init_radii, init_bin_widths = cd.read_distr(outfile, zmax)
    #os.remove(outfile)

    l = sol_str if aerosol==sol else sol_str + " + dust"
    c = "tab:blue" if aerosol==sol else "tab:orange"
    if aerosol == sol:
        ax[0].bar(init_radii, initial_distr, color=c, edgecolor=c, width=init_bin_widths, alpha=0.4, label = l, linewidth=2)
        ax[1].bar(radii, distr1, color=c, edgecolor=c, width=bin_widths, alpha=0.4, label = l, linewidth=2)
    elif aerosol == mix:
        ax[2].bar(init_radii, initial_distr, color=c, edgecolor=c, width=init_bin_widths, alpha=0.4, label = l, linewidth=2)
        ax[3].bar(radii, distr1, color=c, edgecolor=c, width=bin_widths, alpha=0.4, label = l, linewidth=2)


ax[0].set_title('Initial size distribution')
ax[1].set_title('Size distribution at '+str(zmax)+' m')
for axis in ax:
    axis.set_xscale('log')
    axis.set_yscale('log')
ax[0].set_ylabel('droplet concentration [1/mg]')
ax[2].set_ylabel('droplet concentration [1/mg]')               
ax[2].set_xlabel('droplet radius [$\\mu$m]')
ax[3].set_xlabel('droplet radius [$\\mu$m]')

if aerosol_str == "pristine":
    ax[1].set_xlim(10,60)
    ax[3].set_xlim(10,60)
# elif aerosol_str == "polluted":
#     ax[1].set_xlim(6,50)


handles, labels = [], []
for axis in ax:
    h, l = axis.get_legend_handles_labels()
    handles.extend(h)
    labels.extend(l)
by_label = dict(zip(labels, handles))
fig.legend(by_label.values(), by_label.keys(), loc='lower center', bbox_to_anchor=(0.5, 0.89),
           ncol=2, frameon=False)

plt.tight_layout(rect=[0, 0, 1, 0.90])

for axis in [ax[1], ax[3]]:
    xmin, xmax = axis.get_xlim()
    ticks = np.arange(np.ceil(xmin), np.floor(xmax) + 1, 10)
    axis.xaxis.set_major_locator(FixedLocator(ticks))
    axis.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:g}'))
    axis.xaxis.set_minor_locator(FixedLocator([]))
    axis.xaxis.set_minor_formatter(NullFormatter())
    axis.xaxis.get_offset_text().set_visible(False)
    axis.tick_params(axis='x')

plt.savefig(out_png, dpi=200)