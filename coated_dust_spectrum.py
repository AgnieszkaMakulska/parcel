import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")

import numpy as np
from parcel import parcel
from scipy.io import netcdf
import matplotlib.pyplot as plt
from libcloudphxx import common
from matplotlib.ticker import FixedLocator, FuncFormatter, NullFormatter, LogLocator
plt.style.use('seaborn-v0_8')
plt.rcParams.update({
    'font.size': 16,
    'axes.labelsize': 16,
    'axes.titlesize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16
})

from coated_dust import aerosol_spec, run_scheme, read_profiles, read_distr

aerosol_str = "polluted"
rd_insol = 1e-6
rd_sol = 0.02e-6
rd = np.cbrt(rd_insol**3 + rd_sol**3)
epsilon = rd_sol **3 / rd**3

sol = aerosol_spec(aerosol_str, 1.0)
mix = aerosol_spec(aerosol_str, epsilon, rd)


def make_profiles():

    out_png = "plots/rd_insol/single_rd_profile_" + aerosol_str + ".pdf"

    fig, ax = plt.subplots(2, 2, figsize=(8.0, 8.0), sharey=True, squeeze=False)
    ax = ax.flatten()
    
    for aerosol in [sol, mix]:
        outfile = "sol.nc" if aerosol == sol else "mix.nc"
        z, liq_mix_ratio, conc, mean_r, std_dev_area = read_profiles(outfile)

        lw = 2
        if aerosol == sol:
            l = 'sea salt'
            c = 'steelblue'
            s = '-'
        else:
            l = 'sea salt + dust'
            c = 'darkorange'
            s = '-'
        ax[0].plot(liq_mix_ratio, z, color=c, label=l, linestyle=s, linewidth = lw)
        ax[1].plot(conc, z, color=c, label=l, linestyle=s, linewidth = lw)
        ax[2].plot(mean_r, z, color=c, linestyle=s, linewidth = lw)
        ax[3].plot(std_dev_area, z, color=c, label=l, linestyle=s, linewidth = lw)

    ax[0].set_ylabel('z [m]')
    ax[2].set_ylabel('z [m]')
    ax[0].set_xlabel('liquid mix. ratio [g/kg]')
    ax[1].set_xlabel('droplet concentration [1/mg]')
    ax[2].set_xlabel(f'droplet mean radius [$\mu$m]')
    ax[3].set_xlabel(f'std. dev. of droplet area [$\mu$m$^2$]')

    handles, labels = ax[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0.94),
               ncol=2, frameon=False)
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    plt.savefig(out_png, dpi=200)




def make_spectrum():

    out_png = "plots/rd_insol/single_rd_spectrum_" + aerosol_str + ".pdf"
    fig, ax = plt.subplots(1, 2, figsize=(12.0, 5.0), sharey=True, squeeze=False)
    ax = ax.flatten()

    for aerosol in [sol, mix]:
        outfile = "sol.nc" if aerosol == sol else "mix.nc"
        distr1, radii, bin_widths, initial_distr, init_radii, init_bin_widths = read_distr(outfile)

        l = "sea salt" if aerosol==sol else "sea salt + dust"
        c = "tab:blue" if aerosol==sol else "tab:orange"
        ax[0].bar(init_radii, initial_distr, color=c, edgecolor=c, width=init_bin_widths, alpha=0.4, label = l, linewidth=2)
        ax[1].bar(radii, distr1, color=c, edgecolor=c, width=bin_widths, alpha=0.4, label = l, linewidth=2)
    
    ax[0].set_title('Initial size distribution')
    ax[1].set_title('Size distribution at 400 m')
    ax[0].set_xscale('log')
    ax[1].set_xscale('log')
    ax[0].set_yscale('log')
    ax[1].set_yscale('log')
    ax[0].set_ylabel('droplet concentration [1/mg]')        
    ax[0].set_xlabel(f'droplet radius [$\mu$m]')
    ax[1].set_xlabel(f'droplet radius [$\mu$m]')

    if aerosol_str == "pristine":
        ax[1].set_xlim(12.5,17)
    elif aerosol_str == "polluted":
        ax[1].set_xlim(6,12)
    
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


for aerosol in [sol, mix]:
    outfile = "sol.nc" if aerosol == sol else "mix.nc"
    run_scheme(aerosol, outfile)

make_profiles()
make_spectrum()