"""
Checking if mixing between substeps is important
"""

import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")

import numpy as np
from parcel import parcel
from scipy.io import netcdf
import matplotlib.pyplot as plt
from libcloudphxx import common
plt.rcParams.update({'font.size': 16})
from matplotlib.ticker import LogLocator, FuncFormatter, NullFormatter

timesteps = [1]
w_list = [0.25]
z_max_list = [1500]
z1 = 900
z2 = 1200

polluted = '{"polluted": {"kappa": 1.28, "rd_insol" : 0.0, "mean_r": [0.029e-6, 0.071e-6], "gstdev": [1.36, 1.57], "n_tot": [160.0e6, 380.0e6]}}'

pristine = '{"pristine": {"kappa": 1.28, "rd_insol": 0.0, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125.0e6, 65.0e6]}}'

#monomod = '{"monomodal": {"kappa": 1.28, "rd_insol": 0.0, "mean_r": [0.04e-6], "gstdev": [2.2], "n_tot": [1000.0e6]}}'

aerosol_list = [polluted]


def run_scheme(mixing, dt, sstp, aerosol, w_max, z_max, outfile):
    args = dict(
        p_0=90000,
        RH_0=0.97,
        T_0=283,
        aerosol = aerosol,
        #dry_sizes = monodisperse,
        sd_conc=None,
        sd_const_multi=1,
        n_sd_mac=1e20,
        dt=dt,
        z_max=None,
        w = lambda t: w_max if t <= z_max/w_max else -w_max,
        t = z_max / w_max,
        outfile=outfile,
        outfreq=1,
        scheme="lgrngn",
        out_bin='{"liq": {"rght": 1, "moms": [0,1,2,3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0.5e-20},' \
                '"size_distr": {"rght": 20e-6, "moms": [0], "drwt": "wet", "nbin": 100, "lnli": "log", "left": 0.1e-6}}',
        sstp_cond = sstp,
        adaptive_sstp_cond = False,
        sstp_cond_mix   = mixing,
        exact_sstp_cond = True,
        aerosol_independent_of_rhod=True, 
        backend="OpenMP",
        ice_switch = False,
        ice_nucl = False,
        time_dep_ice_nucl = True,
        depo = False
    )
    parcel(**args)


def run(aerosol, w_max, z_max):
    for i in range(len(timesteps)):
        dt = timesteps[i]
        sstp = 10 * dt
        for mixing in [True, False]:
            outfile = "mixing_"+str(mixing)+"_dt_"+str(dt)+".nc"
            run_scheme(mixing, dt, sstp, aerosol, w_max, z_max, outfile)


def read_profiles(outfile):
    with netcdf.netcdf_file(outfile, 'r') as f:
        z = np.array(f.variables['z'][:]).squeeze()
        act_m0 = np.array(f.variables['act_m0'][:]).squeeze()
        act_m1 = np.array(f.variables['act_m1'][:]).squeeze()
        act_m2 = np.array(f.variables['act_m2'][:]).squeeze()
        liq_mix_ratio = np.array(f.variables['liq_m3'][:]).squeeze() * 4/3 * np.pi * common.rho_w
        act_conc = np.array(f.variables['act_m0'][:]).squeeze()
        act_r = np.where(act_conc > 0, np.array(f.variables['act_m1'][:]).squeeze() / np.array(f.variables['act_m0'][:]).squeeze(), 0)
        std_dev_liq = np.sqrt(np.where(act_m0 > 0, 
                        act_m2 / act_m0 - (act_m1 / act_m0)**2, 
                        0))
        rel_disp = np.where(act_r > 0, std_dev_liq / act_r, 0)
    return z/1000, liq_mix_ratio*1e3, act_conc/1e6, act_r*1e6, std_dev_liq*1e6, rel_disp

def read_distr(outfile):
    with netcdf.netcdf_file(outfile, 'r') as f:
        z = np.array(f.variables['z'][:]).squeeze()
        distr = np.array(f.variables['size_distr_m0'][:]).squeeze()
        radii = np.array(f.variables['size_distr_r_wet'][:]).squeeze()
        bin_widths = np.array(f.variables['size_distr_dr_wet'][:]).squeeze()
        distr1 = distr[np.argmin(np.abs(z - z1))]
        distr2 = distr[np.argmin(np.abs(z - z2))]
    return distr1/1e6, distr2/1e6, radii*1e6, bin_widths*1e6

    

def make_figures(aerosol, w_max):

    if aerosol == pristine:
        aerosol_str = "pristine"
    elif aerosol == polluted:
        aerosol_str = "polluted"
    else:
        aerosol_str = "monomodal"
    out_png = "plots/outputs/coupled_uncoupled/"+aerosol_str+"_w_"+str(w_max)

    # plotting profiles
    fig, ax = plt.subplots(len(timesteps), 4, figsize=(16.0, 15.0), sharey=True, squeeze=False)
    for i in range(len(timesteps)):
        dt = timesteps[i]
        for mixing in [True, False]:
            outfile = "mixing_"+str(mixing)+"_dt_"+str(dt)+".nc"
            z, liq_mix_ratio, act_conc, act_r, std_dev_liq, rel_disp = read_profiles(outfile)

            lw = 2
            if mixing:
                l = 'coupled'
                c = 'blue'
                s = '-'
            else:
                l = 'uncoupled'
                c = 'violet'
                s = ':'
            #ax[i,0].plot(liq_mix_ratio, z, color=c, label=l, linestyle=s, linewidth = lw)
            ax[i,0].plot(rel_disp, z, color=c, label=l, linestyle=s, linewidth = lw)
            ax[i,1].plot(act_conc, z, color=c, label=l, linestyle=s, linewidth = lw)
            ax[i,2].plot(act_r, z, color=c, linestyle=s, linewidth = lw)
            ax[i,3].plot(std_dev_liq, z, color=c, label=l, linestyle=s, linewidth = lw)

        ax[i,0].set_ylabel('z [km]')
    ax[-1,0].set_xlabel('rel. disp. of droplet radius')
    #ax[-1,0].set_xlabel('liquid mix. ratio [g/kg]')
    ax[-1,1].set_xlabel('droplet concentration [1/mg]')
    ax[-1,2].set_xlabel(f'droplet mean radius [$\mu$m]')
    ax[-1,3].set_xlabel(f'std. dev. of droplet radius [$\mu$m]')

    handles, labels = ax[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center right", bbox_to_anchor=(0.8, 0.97))
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    dt_labels = ["dt = " + str(dt) + " s" for dt in timesteps]
    for row_idx, label in enumerate(dt_labels):
        axis = ax[row_idx, 0]
        axis.text(0.05, 0.95, label, transform=axis.transAxes, 
                fontsize=18, fontweight='bold', va='top', ha='left',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    plt.suptitle('w = '+str(w_max)+' m/s, '+ aerosol_str)
    plt.savefig(out_png + ".pdf", dpi=200)

    # plotting size distribution
    fig, ax = plt.subplots(len(timesteps), 2, figsize=(16.0, 7.0* len(timesteps)), sharey=True, squeeze=False)
    for i in range(len(timesteps)):
        dt = timesteps[i]
        for mixing in [True, False]:
            outfile = "mixing_"+str(mixing)+"_dt_"+str(dt)+".nc"
            distr1, distr2, radii, bin_widths = read_distr(outfile)
            l = "coupled" if mixing else "uncoupled"
            c = "tab:blue" if mixing else "tab:orange"
            ax[i, 0].bar(radii, distr1, color=c, edgecolor=c, width=bin_widths, alpha=0.4, label = l, linewidth=2)
            ax[i, 1].bar(radii, distr2, color=c, edgecolor=c, width=bin_widths, alpha=0.4, label = l, linewidth=2)
        ax[i,0].set_title(f'z = '+str(z1)+' m')
        ax[i,1].set_title(f'z = '+str(z2)+' m')
        ax[i,0].set_xlim(2,20)
        ax[i,1].set_xlim(2,20)
        ax[i,0].set_xscale('log')
        ax[i,1].set_xscale('log')
        ax[i,0].set_yscale('log')
        ax[i,1].set_yscale('log')
        ax[i, 0].set_ylabel('droplet concentration [1/mg]')        
    ax[-1,0].set_xlabel(f'droplet radius [$\mu$m]')
    ax[-1,1].set_xlabel(f'droplet radius [$\mu$m]')
    ax[0, 0].legend()

    for j in [0, 1]:
        ax[i, j].xaxis.set_major_locator(
            LogLocator(base=10, subs=[1,2,3,4,5,6,7,8,9])
        )
        ax[i, j].xaxis.set_major_formatter(
            FuncFormatter(lambda x, _: f'{x:g}')
        )
        ax[i, j].xaxis.set_minor_formatter(NullFormatter())

    dt_labels = ["dt = " + str(dt) + " s" for dt in timesteps]
    for row_idx, label in enumerate(dt_labels):
        axis = ax[row_idx, 0]
        axis.text(0.05, 0.95, label, transform=axis.transAxes, 
                fontsize=18, fontweight='bold', va='top', ha='left',
                bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
    plt.suptitle('w = '+str(w_max)+' m/s, '+ aerosol_str)
    plt.savefig(out_png + "_distr.pdf", dpi=200)

for w_max,z_max in zip(w_list, z_max_list):
    for aerosol in aerosol_list:
        run(aerosol, w_max, z_max)
        make_figures(aerosol, w_max)
