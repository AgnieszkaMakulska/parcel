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

aerosol_str = "pristine"
zmax = 10

dust = cd.dust_aerosol(epsilon = 0.05)
sol = cd.soluble_aerosol(aerosol_str)


out_png = "plots/rd_insol/dry_spectrum_" + aerosol_str + ".pdf"
fig, ax = plt.subplots(1, 1, figsize=(7.0, 5.0))

for aerosol in [sol, dust]:
    outfile = "dry_spec.nc"
    cd.run_scheme(aerosol, zmax, outfile, outfreq = 10, spec = True)
    distr, radii, bin_widths = cd.read_dry_distr(outfile)
    os.remove(outfile)

    l = "ammonium sulfate" if aerosol==sol else "coated dust"
    ax.bar(radii, distr, width=bin_widths, alpha=0.5, label = l, linewidth=2)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylabel('concentration [1/mg]')        
ax.set_xlabel('dry radius [$\\mu$m]')


handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0.89),
            ncol=2, frameon=False)
plt.tight_layout(rect=[0, 0, 1, 0.90])
plt.savefig(out_png, dpi=200)