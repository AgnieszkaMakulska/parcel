"""
Testing adaptive timesteps with the WBF process
"""

import sys, os
sys.path.insert(0, "../")
sys.path.insert(0, "./")

import numpy as np
from parcel import parcel
from scipy.io import netcdf
import matplotlib.pyplot as plt
from functions import rh_to_rh_i
from libcloudphxx import common
from matplotlib.ticker import MaxNLocator
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection

sstp_max = 10
w_list = [2.5]#[1., 2.5, 5.]
z_max_list = [2000]#[1250., 2000, 3000.]

polluted = '{"polluted": {"kappa": 0.61, "rd_insol" : 0.0, "mean_r": [0.029e-6, 0.071e-6], "gstdev": [1.36, 1.57], "n_tot": [160.0e6, 380.0e6]},' \
                '"INP": {"kappa": 0.61, "rd_insol" : 0.5e-6, "mean_r": [0.029e-6], "gstdev": [1.36], "n_tot": [10.0e6]}}' # low concentration of INPs

pristine = '{"pristine": {"kappa": 0.61, "rd_insol": 0.0, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125.0e6, 65.0e6]},' \
                '"INP": {"kappa": 0.61, "rd_insol" : 0.5e-6, "mean_r": [0.029e-6], "gstdev": [1.36], "n_tot": [10.0e6]}}' # low concentration of INPs

def run_scheme(outfile, aerosol, w_max, z_max, adaptive, epsilon):
    args = dict(
        p_0=100000,
        RH_0=0.9,
        T_0=277,
        aerosol = aerosol,
        sd_conc=500,
        dt=1,
        z_max=None,
        #w=lambda t: w_max * np.pi / 2. * np.sin(np.pi*t*w_max/z_max),
        w = lambda t: w_max if t <= z_max/w_max else -w_max,
        t = 2 * z_max / w_max,
        outfile=outfile,
        outfreq=10,
        scheme="lgrngn",
        out_bin='{"liq": {"rght": 1, "moms": [0,1,2,3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0.5e-20},' \
                '"ice": {"rght": 1, "moms": [0,1,2,3], "drwt": "ice_a", "nbin": 1, "lnli": "lin", "left": 0.5e-20}}',
        sstp_cond = sstp_max,
        adaptive_sstp_cond = adaptive,
        sstp_cond_mix   = False,
        exact_sstp_cond = True,
        aerosol_independent_of_rhod=True, 
        backend="OpenMP",
        ice_switch = True,
        ice_nucl = True,
        time_dep_ice_nucl = False,
        depo = True,
        sstp_cond_adapt_drw2_eps = epsilon
        #sstp_cond_adapt_drw2_max=100,
        #sstp_cond_act=1
    )

    parcel(**args)

    with netcdf.netcdf_file(outfile, 'r') as f:
        z = np.array(f.variables['z'][:]).squeeze()
        RH = np.array(f.variables['RH'][:]).squeeze()
        T = np.array(f.variables['T'][:]).squeeze()
        rv = np.array(f.variables['r_v'][:]).squeeze()
        liq_m0 = np.array(f.variables['liq_m0'][:]).squeeze()
        liq_m1 = np.array(f.variables['liq_m1'][:]).squeeze()
        liq_m2 = np.array(f.variables['liq_m2'][:]).squeeze()
        ice_m0 = np.array(f.variables['ice_m0'][:]).squeeze()
        ice_m1 = np.array(f.variables['ice_m1'][:]).squeeze()
        ice_m2 = np.array(f.variables['ice_m2'][:]).squeeze()

        ice_mix_ratio = np.array(f.variables['ice_mix_ratio'][:]).squeeze()
        liq_mix_ratio = np.array(f.variables['liq_m3'][:]).squeeze() * 4/3 * np.pi * common.rho_w  #multiply by density of water
        ice_conc = np.array(f.variables['ice_m0'][:]).squeeze()  # 1/kg
        act_conc = np.array(f.variables['act_m0'][:]).squeeze()
        ice_r = np.where(ice_conc > 0, np.array(f.variables['ice_m1'][:]).squeeze() / np.array(f.variables['ice_m0'][:]).squeeze(), 0)
        act_r = np.where(act_conc > 0, np.array(f.variables['act_m1'][:]).squeeze() / np.array(f.variables['act_m0'][:]).squeeze(), 0)
        sstp_cond_mean = np.array(f.variables['sstp_cond_mean'][:]) if 'sstp_cond_mean' in f.variables else None

        variance_liq = np.where(liq_m0 > 0, 
                           liq_m2 / liq_m0 - (liq_m1 / liq_m0)**2, 
                           0)
        variance_ice = np.where(ice_m0 > 0, 
                           ice_m2 / ice_m0 - (ice_m1 / ice_m0)**2, 
                           0)

        if sstp_cond_mean is not None:
            sstp_cond_mean[0] = sstp_cond_mean[1] # at t=0 sstp_cond_mean=0, because its set only during the firs step (?)
        sstp_dep_mean = np.array(f.variables['sstp_dep_mean'][:]) if 'sstp_dep_mean' in f.variables else None
        if sstp_dep_mean is not None:
            sstp_dep_mean[0] = sstp_dep_mean[1] # at t=0 sstp_cond_mean=0, because its set only during the firs step (?)
    return RH, T, rv, z, ice_mix_ratio, liq_mix_ratio, ice_conc, act_conc, ice_r, act_r, sstp_cond_mean, sstp_dep_mean, variance_liq, variance_ice



def make_figure(aerosol, w_max, z_max):

    fig, ax = plt.subplots(2, 5, figsize=(12.0, 12.0), sharey=True, squeeze=True)

    #for adaptive, epsilon in [(True, 1e-1),(True, 1e-2), (True, 1e-3), (False, None)]:
    for adaptive, epsilon in [(True, 1e-1)]:

        outfile = f"test_WBF.nc"
        RH, T, rv, z, ice_mix_ratio, liq_mix_ratio, ice_conc, act_conc, ice_r, act_r, sstp_cond_mean, sstp_dep_mean, variance_liq, variance_ice = run_scheme(outfile, aerosol, w_max, z_max, adaptive, epsilon)           

        if not adaptive:
            c = 'lightgrey'
            s = ':'
            l = 'non-adaptive'
        else:
            s = '-'
            l = '$\epsilon$ = '+str(epsilon)
            if epsilon == 1e-3:
                c = 'blue'
            elif epsilon == 1e-2:
                c = 'purple'
            else:
                c = 'violet'

        m = len(RH)//2

        ax[1,0].plot(ice_mix_ratio * 1e3, z, color=c, label=l, linestyle=s)
        ax[0,0].plot(liq_mix_ratio * 1e3, z, color=c, label=l, linestyle=s)
        ax[1,1].plot(ice_conc / 1e6, z, color=c, label=l, linestyle=s)
        ax[0,1].plot(act_conc / 1e6, z, color=c, label=l, linestyle=s)
        ax[1,2].plot(ice_r * 1e6, z, color=c, linestyle=s)
        ax[0,2].plot(act_r * 1e6, z, color=c, linestyle=s)
        ax[1,3].plot(np.sqrt(variance_ice)*1e6, z, color=c, label=l, linestyle=s)
        ax[0,3].plot(np.sqrt(variance_liq)*1e6, z, color=c, label=l, linestyle=s)
        if adaptive:
            ax[1,4].plot(sstp_dep_mean, z, color=c, linestyle=s)
            ax[0,4].plot(sstp_cond_mean, z, color=c, linestyle=s)


        x_start = (liq_mix_ratio[:m] * 1e3)[len(liq_mix_ratio[:m])//2]
        y_start = z[:m][len(z[:m])//2]
        x_end = (liq_mix_ratio[:m] * 1e3)[len(liq_mix_ratio[:m])//2+2]
        y_end = z[:m][len(z[:m])//2+2]
        ax[0,0].annotate(
            '',
            xy=(x_end, y_end),    # Grot strzałki
            xytext=(x_start, y_start), # Początek strzałki
            arrowprops=dict(
                arrowstyle='->',
                color=c,
                linewidth=2
            )
        )


        x_end = (liq_mix_ratio[m:] * 1e3)[len(liq_mix_ratio[m:])//2]
        y_end = z[:m][len(z[m:])//2]
        x_start= (liq_mix_ratio[m:] * 1e3)[len(liq_mix_ratio[m:])//2+2]
        y_start = z[:m][len(z[m:])//2+2]
        ax[0,0].annotate(
            '',
            xy=(x_end, y_end),    # Grot strzałki
            xytext=(x_start, y_start), # Początek strzałki
            arrowprops=dict(
                arrowstyle='->',
                color=c,
                linewidth=2
            )
        )

    ax[0,3].legend()
    handles, labels = ax[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.15, 0.2))

    ax[0,0].set_ylabel('z [m]')
    ax[1,0].set_ylabel('z [m]')
    ax[0,0].set_xlabel('liq mix ratio [g/kg]')
    ax[1,0].set_xlabel('ice mix ratio [g/kg]')
    ax[0,1].set_xlabel('liq concentration [1/mg]')
    ax[1,1].set_xlabel('ice concentration [1/mg]')
    ax[0,2].set_xlabel("liq average radius [um]")
    ax[1,2].set_xlabel("ice average radius [um]")
    ax[0,3].set_xlabel('liq standard deviation [um]')
    ax[1,3].set_xlabel('ice standard deviation [um]')
    ax[0,4].set_xlabel("liq sstp mean")
    ax[1,4].set_xlabel("ice sstp mean")

    # ax[0,3].ticklabel_format(style='plain', useOffset=False)
    # ax[0,3].xaxis.set_major_locator(MaxNLocator(4))

    ax[0,4].set_xlim(0,10.1)
    ax[1,4].set_xlim(0,10.1)

    fig.tight_layout(rect=(0, 0.10, 1, 0.97))
    aerosol_str = "pristine" if aerosol==pristine else "polluted"
    out_png = "test_adaptive_WBF_w_"+str(w_max)+"_"+aerosol_str+".png"
    plt.suptitle('w = '+str(w_max)+' m/s, '+aerosol_str)
    plt.savefig(out_png, dpi=200)

    return fig

for w_max,z_max in zip(w_list, z_max_list):
    print(w_max, z_max)
    for aerosol in [pristine]: #, polluted]:
        make_figure(aerosol, w_max, z_max)
plt.show()