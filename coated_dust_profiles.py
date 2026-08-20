import coated_dust as cd
import sys
import os
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8')
plt.rcParams.update({
    'font.size': 16,
    'axes.labelsize': 16,
    'axes.titlesize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16,
    'lines.linewidth': 1.5
})


aerosol_str = "pristine"
out_png = "plots/rd_insol/profiles_" + aerosol_str + ".pdf"

rd = 3e-6
epsilon_list = [0.001, 0.05, 0.1, 0.9]

mosaic = [
    ["ax1", "ax1", "ax2", "ax2"],
    ["ax3", "ax3", "ax4", "ax4"],
    [".", "ax5", "ax5", "."],
]
fig, axd = plt.subplot_mosaic(mosaic, figsize=(8.0, 11.0), sharey=True)
ax = list(axd.values())

# no dust
sol = cd.soluble_aerosol(aerosol_str)
l = 'no dust'
outfile = "sol.nc"
cd.run_scheme(sol, outfile, outfreq = 5)
z, rh, liq_mix_ratio, conc, mean_r, std_dev_r = cd.read_profiles(outfile)
os.remove(outfile)
ax[0].plot(liq_mix_ratio, z, label=l, color="black",linestyle="--")
ax[1].plot(conc, z, label=l, color="black",linestyle="--")
ax[2].plot(mean_r, z, label=l, color="black",linestyle="--")
ax[3].plot(std_dev_r, z, label=l, color="black",linestyle="--")
ax[4].plot((rh-1)*100, z, label=l, color="black",linestyle="--")

# coated dust
for epsilon in epsilon_list:
    aerosol = cd.mixed_aerosol(aerosol_str, epsilon, rd)
    l = '$\\epsilon$ = ' + str(round(epsilon, 5))
    outfile = str(epsilon)+".nc"
    cd.run_scheme(aerosol, outfile, outfreq = 5)
    z, rh, liq_mix_ratio, conc, mean_r, std_dev_r = cd.read_profiles(outfile)
    os.remove(outfile)
    ax[0].plot(liq_mix_ratio, z, label=l)
    ax[1].plot(conc, z, label=l)
    ax[2].plot(mean_r, z, label=l)
    ax[3].plot(std_dev_r, z, label=l)
    ax[4].plot((rh-1)*100, z, label=l)

ax[0].set_ylabel('z [m]')
ax[2].set_ylabel('z [m]')
ax[4].set_ylabel('z [m]')
ax[0].set_xlabel('liquid mix. ratio [g/kg]')
ax[1].set_xlabel('droplet concentration [1/mg]')
ax[2].set_xlabel('droplet mean radius [$\\mu$m]')
ax[3].set_xlabel('std. dev. of droplet radius [$\\mu$m]')
ax[4].set_xlabel('RH [%]')

if aerosol_str == "pristine":
    ax[1].set_xlim(57,61)
    ax[4].set_xlim(0.0,0.8)
#     ax[2].set_xlim(2,10)
# elif aerosol_str == "polluted":
#     ax[1].set_xlim(230,370)
#     ax[2].set_xlim(2,9)


handles, labels = ax[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0.89),
            ncol=2, frameon=False)
ax[4].tick_params(labelleft=True)
plt.tight_layout(rect=[0, 0, 1, 0.9])
plt.savefig(out_png, dpi=200, bbox_inches='tight')