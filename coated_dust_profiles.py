import coated_dust as cd
import sys
import os
import matplotlib.pyplot as plt
plt.style.use('seaborn-v0_8')
plt.rcParams.update({
    'font.size': 14,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 14,
    'lines.linewidth': 1.5
})


aerosol_str = "pristine"
out_png = "plots/rd_insol/profiles_" + aerosol_str + ".pdf"
zmax = 400

epsilon_list = [0.001, 0.05, 0.5, 1.]

mosaic = [
    ["ax1", "ax1", "ax2", "ax2"],
    ["ax3", "ax3", "ax4", "ax4"],
    [".", "ax5", "ax5", "."],
]
fig, axd = plt.subplot_mosaic(mosaic, figsize=(8.0, 11.0), sharey=False)
ax = list(axd.values())

# no dust
sol = cd.soluble_aerosol(aerosol_str)
l = 'sulfate'
outfile = "sol.nc"
cd.run_scheme(sol, zmax, outfile, outfreq = 5)
z, rh, liq_mix_ratio, conc, mean_r, std_dev_r = cd.read_profiles(outfile)
os.remove(outfile)
ax[0].plot(conc, z, label=l, color="black",linestyle="--")
ax[1].plot(liq_mix_ratio, z, label=l, color="black",linestyle="--")
ax[2].plot(mean_r, z, label=l, color="black",linestyle="--")
ax[3].plot(std_dev_r, z, label=l, color="black",linestyle="--")
ax[4].plot((rh-1)*100, z, label=l, color="black",linestyle="--")

# coated dust
for epsilon in epsilon_list:
    aerosol = cd.mixed_aerosol(aerosol_str, epsilon)
    if epsilon == 1.0:
        l = 'sea salt ($\\epsilon$ = 1) + sulfate'
    else:
        l = 'dust ($\\epsilon$ = ' + str(round(epsilon, 5)) + ') + sulfate'
    outfile = str(epsilon)+".nc"
    cd.run_scheme(aerosol, zmax, outfile, outfreq = 1)
    z, rh, liq_mix_ratio, conc, mean_r, std_dev_r = cd.read_profiles(outfile)
    #os.remove(outfile)
    ax[0].plot(conc, z, label=l)
    ax[1].plot(liq_mix_ratio, z, label=l)
    ax[2].plot(mean_r, z, label=l)
    ax[3].plot(std_dev_r, z, label=l)
    ax[4].plot((rh-1)*100, z, label=l)

ax[0].set_ylabel('z [m]')
ax[2].set_ylabel('z [m]')
ax[4].set_ylabel('z [m]')
ax[0].set_xlabel('droplet concentration [1/mg]')
ax[1].set_xlabel('liquid mix. ratio [g/kg]')
ax[2].set_xlabel('droplet mean radius [$\\mu$m]')
ax[3].set_xlabel('std. dev. of droplet radius [$\\mu$m]')
ax[4].set_xlabel('RH [%]')

if aerosol_str == "pristine":
    ax[0].set_xlim(57.5,63.) # w 1
    ax[4].set_xlim(0.6,0.8)
    ax[4].set_ylim(50,120)
    #ax[0].set_xlim(57,100) # w 2.5
    #ax[1].set_xlim(0.0,1.5)
    # ax[0].set_xlim(50,60) # w 0.5
    # ax[1].set_xlim(0.0,1.)

elif aerosol_str == "polluted":
    ax[0].set_xlim(390,400) # w 1
    ax[4].set_xlim(0.2,0.4)
    ax[4].set_ylim(50,100)


handles, labels = ax[0].get_legend_handles_labels()
ax[4].tick_params(labelleft=True)
# fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0.89),
#             ncol=2, frameon=False)
fig.legend(handles[:1], labels[:1], loc='lower center', bbox_to_anchor=(0.5, 0.90),
           frameon=False)
fig.legend(handles[1:], labels[1:], loc='lower center', bbox_to_anchor=(0.5, 0.83),
           ncol=2, frameon=False)
#plt.tight_layout(rect=[0, 0, 1, 0.9])
plt.tight_layout(rect=[0, 0, 1, 0.84])
plt.savefig(out_png, dpi=200, bbox_inches='tight')