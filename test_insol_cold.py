"""
Checking if insoluble component is important for condensation and deposition
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

w_list = [1.]
z_max_list = [4000]
z1 = 1000

#composition of the mixed particle
rd_insol = 0.5e-6
rd_sol = 0.02e-6
rd = np.cbrt(rd_insol**3 + rd_sol**3)
epsilon = rd_sol **3 / rd**3
print(epsilon)

sol = f'{{"soluble": {{"kappa": 1.28, "sol_frac": 1.0, "mean_r": [{rd_sol}], "gstdev": [1.4], "n_tot": [200.0e6]}} }}'

mix = f'{{"soluble": {{"kappa": 1.28, "sol_frac": 1.0, "mean_r": [{rd_sol}], "gstdev": [1.4], "n_tot": [199.0e6]}}, \
           "mixed": {{"kappa": 1.28, "sol_frac": {epsilon}, "mean_r": [{rd}], "gstdev": [1.4], "n_tot": [1.0e6]}} }}'


def run_scheme(aerosol, w, z_max, outfile):
    args = dict(
        p_0=90000,
        RH_0=0.97,
        T_0=260,
        aerosol = aerosol,
        w = w,
        #dry_sizes = monodisperse,
        sd_conc=100,
        #sd_const_multi=1000000,
        #n_sd_max=1e7,
        dt=1,
        z_max = z_max,
        outfile=outfile,
        outfreq=1,
        scheme="lgrngn",
        out_bin='{"liq": {"rght": 1, "moms": [0,1,2,3,4], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0.5e-20},' \
                '"ice": {"rght": 1, "moms": [0,1,2,3,4], "drwt": "ice_a", "nbin": 1, "lnli": "lin", "left": 0.5e-20},' \
                '"initial_spec": {"rght": 3e-6, "moms": [0], "drwt": "wet", "nbin": 100, "lnli": "log", "left": 0.01e-6},' \
                '"spec": {"rght": 30e-6, "moms": [0], "drwt": "wet", "nbin": 1000, "lnli": "log", "left": 1e-6}}',

        sstp_cond = 10,
        ice_switch = True,
        ice_nucl = True,
        time_dep_ice_nucl = True,
        depo = True
    )
    parcel(**args)


def read_profiles(outfile):
    with netcdf.netcdf_file(outfile, 'r') as f:
        z = np.array(f.variables['z'][:]).squeeze()
        T = np.array(f.variables['T'][:]).squeeze()
        act_m0 = np.array(f.variables['act_m0'][:]).squeeze()
        act_m1 = np.array(f.variables['act_m1'][:]).squeeze()
        act_m2 = np.array(f.variables['act_m2'][:]).squeeze()
        act_m4 = np.array(f.variables['act_m4'][:]).squeeze()
        liq_m3 = np.array(f.variables['liq_m3'][:]).squeeze()
        ice_m0 = np.array(f.variables['ice_m0'][:]).squeeze()
        ice_m1 = np.array(f.variables['ice_m1'][:]).squeeze()
        ice_m2 = np.array(f.variables['ice_m2'][:]).squeeze()
        ice_m3 = np.array(f.variables['ice_m3'][:]).squeeze()
        ice_m4 = np.array(f.variables['ice_m4'][:]).squeeze()
        liq_mix_ratio = liq_m3 * 4/3 * np.pi * common.rho_w
        conc = act_m0
        mean_r = np.where(act_m0 > 0, act_m1 / act_m0, 0)
        std_dev_area = np.sqrt(np.where(act_m0 > 0, 
                        act_m4 / act_m0 - (act_m2 / act_m0)**2, 
                        0))
        ice_mix_ratio = ice_m3 * 4/3 * np.pi * common.rho_i
        ice_conc = ice_m0
        ice_mean_r = np.where(ice_m0 > 0, ice_m1 / ice_m0, 0)
        ice_std_dev_area = np.sqrt(np.where(ice_m0 > 0, 
                        ice_m4 / ice_m0 - (ice_m2 / ice_m0)**2, 
                        0))
    return z/1000, liq_mix_ratio*1e3, conc/1e6, mean_r*1e6, std_dev_area*1e12, ice_mix_ratio*1e3, ice_conc/1e6, ice_mean_r*1e6, ice_std_dev_area*1e12, T

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
        distr1 = distr[np.argmin(np.abs(z - z1))]
        #distr2 = distr[np.argmin(np.abs(z - z2))]
    return distr1/1e6, radii*1e6, bin_widths*1e6, initial_distr/1e6, init_radii*1e6, init_bin_widths*1e6


def make_profiles(w):

    out_png = "plots/rd_insol/w_"+str(w)

    fig, ax = plt.subplots(2, 4, figsize=(16.0, 15.0), sharey=True, squeeze=False)
    
    for aerosol in [sol, mix]:
        if aerosol == sol:
            aerosol_str = "sol"
        elif aerosol == mix:
            aerosol_str = "mix"
        outfile = aerosol_str+".nc"
        z, liq_mix_ratio, conc, mean_r, std_dev_area, ice_mix_ratio, ice_conc, ice_mean_r, ice_std_dev_area, T = read_profiles(outfile)

        lw = 2
        if aerosol == sol:
            l = 'soluble'
            c = 'blue'
            s = '-'
        else:
            l = 'mixed'
            c = 'violet'
            s = ':'
        ax[0,0].plot(liq_mix_ratio, z, color=c, label=l, linestyle=s, linewidth = lw)
        ax[0,1].plot(conc, z, color=c, label=l, linestyle=s, linewidth = lw)
        ax[0,2].plot(mean_r, z, color=c, linestyle=s, linewidth = lw)
        ax[0,3].plot(std_dev_area, z, color=c, label=l, linestyle=s, linewidth = lw)
        ax[1,0].plot(ice_mix_ratio, z, color=c, label=l, linestyle=s, linewidth = lw)
        ax[1,1].plot(ice_conc, z, color=c, label=l, linestyle=s, linewidth = lw)
        ax[1,2].plot(ice_mean_r, z, color=c, linestyle=s, linewidth = lw)
        #ax[1,3].plot(ice_std_dev_area, z, color=c, label=l, linestyle=s, linewidth = lw)
        ax[1,3].plot(T, z, color=c, label=l, linestyle=s, linewidth = lw)

    ax[0,0].set_ylabel('z [km]')
    ax[0,0].set_xlabel('liquid mix. ratio [g/kg]')
    ax[0,1].set_xlabel('droplet concentration [1/mg]')
    ax[0,2].set_xlabel(f'droplet mean radius [$\mu$m]')
    ax[0,3].set_xlabel(f'std. dev. of droplet area [$\mu$m$^2$]')
    ax[1,0].set_xlabel('ice mix. ratio [g/kg]')
    ax[1,1].set_xlabel('ice concentration [1/mg]')
    ax[1,2].set_xlabel(f'ice mean radius [$\mu$m]')
    #ax[1,3].set_xlabel(f'std. dev. of ice area [$\mu$m$^2$]')
    ax[1,3].set_xlabel('T [K]')

    handles, labels = ax[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center right", bbox_to_anchor=(0.8, 0.97))
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    plt.suptitle('w = '+str(w)+' m/s')
    plt.savefig(out_png + ".pdf", dpi=200)




def make_spectrum(w):

    out_png = "plots/rd_insol/w_"+str(w)
    fig, ax = plt.subplots(1, 2, figsize=(16.0, 7.0), sharey=True, squeeze=False)

    for aerosol in [sol, mix]:
        if aerosol == sol:
            aerosol_str = "sol"
        elif aerosol == mix:
            aerosol_str = "mix"
        outfile = aerosol_str+".nc"
        distr1, radii, bin_widths, initial_distr, init_radii, init_bin_widths = read_distr(outfile)

        l = "soluble" if aerosol==sol else "mixed"
        c = "tab:blue" if aerosol==sol else "tab:orange"
        ax[0, 0].bar(init_radii, initial_distr, color=c, edgecolor=c, width=init_bin_widths, alpha=0.4, label = l, linewidth=2)
        ax[0, 1].bar(radii, distr1, color=c, edgecolor=c, width=bin_widths, alpha=0.4, label = l, linewidth=2)
    
    ax[0,0].set_title(f'z = 0 m')
    ax[0,1].set_title(f'z = '+str(z1)+' m')
    # ax[0,0].set_xlim(0,5)
    ax[0,1].set_xlim(14,20)
    ax[0,0].set_xscale('log')
    ax[0,1].set_xscale('log')
    ax[0,0].set_yscale('log')
    ax[0,1].set_yscale('log')
    ax[0, 0].set_ylabel('droplet concentration [1/mg]')        
    ax[0,0].set_xlabel(f'droplet radius [$\mu$m]')
    ax[0,1].set_xlabel(f'droplet radius [$\mu$m]')
    ax[0, 1].legend()

    # for j in [0, 1]:
    #     ax[0, j].xaxis.set_major_locator(
    #         LogLocator(base=10, subs=[1,2,3,4,5,6,7,8,9])
    #     )
    #     ax[0, j].xaxis.set_major_formatter(
    #         FuncFormatter(lambda x, _: f'{x:g}')
    #     )
    #     ax[0, j].xaxis.set_minor_formatter(NullFormatter())
    plt.suptitle('w = '+str(w)+' m/s')
    plt.savefig(out_png + "_distr.pdf", dpi=200)


for w,z_max in zip(w_list, z_max_list):
    for aerosol in [sol, mix]:
        if aerosol == sol:
            aerosol_str = "sol"
        elif aerosol == mix:
            aerosol_str = "mix"
        outfile = aerosol_str+".nc"
        run_scheme(aerosol, w, z_max, outfile)
    make_profiles(w)
    #make_spectrum(w)