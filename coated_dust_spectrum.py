"""
Checking if insoluble component is important for condensation
"""

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
plt.rcParams.update({'font.size': 16})
plt.rcParams.update({
    'axes.labelsize': 16,
    'axes.titlesize': 16,
    'xtick.labelsize': 16,
    'ytick.labelsize': 16,
    'legend.fontsize': 16
})

w = 1.
z_max = 400
z_distr = 400

#composition of the mixed particle
rd_insol = 1e-6
rd_sol = 0.02e-6
rd = np.cbrt(rd_insol**3 + rd_sol**3)
epsilon = rd_sol **3 / rd**3

# sol = '{"polluted": {"kappa": 1.28, "sol_frac": 1.0, "mean_r": [0.029e-6, 0.071e-6], "gstdev": [1.36, 1.57], "n_tot": [160.0e6, 380.0e6]}}'
# mix = f'{{"polluted": {{"kappa": 1.28, "sol_frac": 1.0, "mean_r": [0.029e-6, 0.071e-6], "gstdev": [1.36, 1.57], "n_tot": [160.0e6, 380.0e6]}}, \
#            "mixed": {{"kappa": 1.28, "sol_frac": {epsilon}, "mean_r": [{rd}], "gstdev": [1.4], "n_tot": [1.0e6]}} }}'


sol = '{"pristine": {"kappa": 1.28, "sol_frac": 1.0, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125.0e6, 65.0e6]}}'
mix = f'{{"pristine": {{"kappa": 1.28, "sol_frac": 1.0, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125.0e6, 65.0e6]}}, \
           "mixed": {{"kappa": 1.28, "sol_frac": {epsilon}, "mean_r": [{rd}], "gstdev": [1.4], "n_tot": [1.0e6]}} }}'

def run_scheme(aerosol, outfile):
    args = dict(
        p_0=90000,
        RH_0=0.97,
        T_0=283,
        aerosol = aerosol,
        w = w,
        sd_conc=1000,
        #sd_const_multi=1000000,
        #n_sd_max=1e7,
        dt = 1,
        z_max = z_max,
        outfile=outfile,
        outfreq=1,
        scheme="lgrngn",
        out_bin='{"liq": {"rght": 1, "moms": [0,1,2,3,4], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0.5e-20},' \
                '"initial_spec": {"rght": 3e-6, "moms": [0], "drwt": "wet", "nbin": 100, "lnli": "log", "left": 0.01e-6},' \
                '"spec": {"rght": 30e-6, "moms": [0], "drwt": "wet", "nbin": 1000, "lnli": "log", "left": 1e-6}}',

        sstp_cond = 10,
        ice_switch = False,
        ice_nucl = False,
        time_dep_ice_nucl = True,
        depo = False,
        backend = "gpu"
    )
    parcel(**args)


def read_profiles(outfile):
    with netcdf.netcdf_file(outfile, 'r') as f:
        z = np.array(f.variables['z'][:]).squeeze()
        act_m0 = np.array(f.variables['act_m0'][:]).squeeze()
        act_m1 = np.array(f.variables['act_m1'][:]).squeeze()
        act_m2 = np.array(f.variables['act_m2'][:]).squeeze()
        act_m4 = np.array(f.variables['act_m4'][:]).squeeze()
        liq_m0 = np.array(f.variables['liq_m0'][:]).squeeze()
        liq_m1 = np.array(f.variables['liq_m1'][:]).squeeze()
        liq_m2 = np.array(f.variables['liq_m2'][:]).squeeze()
        liq_m3 = np.array(f.variables['liq_m3'][:]).squeeze()
        liq_m4 = np.array(f.variables['liq_m4'][:]).squeeze()
        liq_mix_ratio = liq_m3 * 4/3 * np.pi * common.rho_w
        conc = act_m0
        mean_r = np.where(act_m0 > 0, act_m1 / act_m0, 0)
        std_dev_area = np.sqrt(np.where(act_m0 > 0, 
                        act_m4 / act_m0 - (act_m2 / act_m0)**2, 
                        0)) * 4 * np.pi
    return z, liq_mix_ratio*1e3, conc/1e6, mean_r*1e6, std_dev_area*1e12

def read_distr(outfile):
    with netcdf.netcdf_file(outfile, 'r') as f:
        z = np.array(f.variables['z'][:]).squeeze()
        distr = np.array(f.variables['spec_m0'][:]).squeeze()
        radii = np.array(f.variables['spec_r_wet'][:]).squeeze()
        bin_widths = np.array(f.variables['spec_dr_wet'][:]).squeeze()
        init_distr = np.array(f.variables['initial_spec_m0'][:]).squeeze()
        init_radii = np.array(f.variables['initial_spec_r_wet'][:]).squeeze()
        init_bin_widths = np.array(f.variables['initial_spec_dr_wet'][:]).squeeze()
        initial_distr = init_distr[np.argmin(np.abs(z))]
        distr1 = distr[np.argmin(np.abs(z - z_distr))]
    return distr1/1e6, radii*1e6, bin_widths*1e6, initial_distr/1e6, init_radii*1e6, init_bin_widths*1e6


def make_profiles():

    out_png = "plots/rd_insol/single_rd_profile.pdf"

    fig, ax = plt.subplots(2, 2, figsize=(8.0, 8.0), sharey=True, squeeze=False)
    ax = ax.flatten()
    
    for aerosol in [sol, mix]:
        if aerosol == sol:
            aerosol_str = "sol"
        elif aerosol == mix:
            aerosol_str = "mix"
        outfile = aerosol_str+".nc"
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

    out_png = "plots/rd_insol/single_rd_spectrum.pdf"
    fig, ax = plt.subplots(1, 2, figsize=(12.0, 5.0), sharey=True, squeeze=False)
    ax = ax.flatten()

    for aerosol in [sol, mix]:
        if aerosol == sol:
            aerosol_str = "sol"
        elif aerosol == mix:
            aerosol_str = "mix"
        outfile = aerosol_str+".nc"
        distr1, radii, bin_widths, initial_distr, init_radii, init_bin_widths = read_distr(outfile)

        l = "sea salt" if aerosol==sol else "sea salt + dust"
        c = "tab:blue" if aerosol==sol else "tab:orange"
        ax[0].bar(init_radii, initial_distr, color=c, edgecolor=c, width=init_bin_widths, alpha=0.4, label = l, linewidth=2)
        ax[1].bar(radii, distr1, color=c, edgecolor=c, width=bin_widths, alpha=0.4, label = l, linewidth=2)
    
    ax[0].set_title('Initial size distribution')
    ax[1].set_title('Size distribution at '+str(z_distr)+' m')
    ax[1].set_xlim(12.5,17)
    ax[0].set_xscale('log')
    ax[1].set_xscale('log')
    ax[0].set_yscale('log')
    ax[1].set_yscale('log')
    ax[0].set_ylabel('droplet concentration [1/mg]')        
    ax[0].set_xlabel(f'droplet radius [$\mu$m]')
    ax[1].set_xlabel(f'droplet radius [$\mu$m]')
    
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
    if aerosol == sol:
        aerosol_str = "sol"
    elif aerosol == mix:
        aerosol_str = "mix"
    outfile = aerosol_str + ".nc"
    run_scheme(aerosol, outfile)
    make_profiles()
    make_spectrum()