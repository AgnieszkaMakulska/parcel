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
from libcloudphxx import common
plt.rcParams.update({'font.size': 14})

w_list = [2.5]#[1., 2.5, 5.]
z_max_list = [4000]#[1250., 2500, 3000.]
outfile = f"test_WBF.nc"

polluted = '{"polluted": {"kappa": 0.61, "rd_insol" : 0.0, "mean_r": [0.029e-6, 0.071e-6], "gstdev": [1.36, 1.57], "n_tot": [160.0e6, 380.0e6]},' \
                '"INP": {"kappa": 0.61, "rd_insol" : 0.5e-6, "mean_r": [0.029e-6], "gstdev": [1.36], "n_tot": [10.0e6]}}' # low concentration of INPs

pristine = '{"pristine": {"kappa": 0.61, "rd_insol": 0.5e-6, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125.0e6, 65.0e6]},' \
                '"INP": {"kappa": 0.61, "rd_insol" : 0.5e-6, "mean_r": [0.029e-6], "gstdev": [1.36], "n_tot": [10.0e6]}}' # low concentration of INPs

aerosol_list = [polluted]

def run_scheme(aerosol, w_max, z_max, adaptive, epsilon, sstp_act, sstp_cond):
    args = dict(
        p_0=100000,
        RH_0=0.8,
        T_0=265.,
        aerosol = aerosol,
        sd_conc=100,
        dt=1,
        z_max=None,
        w = lambda t: w_max if t <= z_max/w_max else -w_max,
        t = 2 * z_max / w_max,
        outfile=outfile,
        outfreq=10,
        scheme="lgrngn",
        out_bin='{"liq": {"rght": 1, "moms": [0,1,2,3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0.5e-20},' \
                '"ice": {"rght": 1, "moms": [0,1,2,3], "drwt": "ice_a", "nbin": 1, "lnli": "lin", "left": 0.5e-20}}',
        sstp_cond = sstp_cond,
        adaptive_sstp_cond = adaptive,
        sstp_cond_mix   = False,
        exact_sstp_cond = True,
        aerosol_independent_of_rhod=True, 
        backend="OpenMP",
        ice_switch = True,
        ice_nucl = True,
        time_dep_ice_nucl = True,
        depo = True,
        sstp_cond_adapt_drw2_eps = epsilon,
        sstp_cond_act = sstp_act
        #sstp_cond_adapt_drw2_max=100,
    )

    parcel(**args)

    with netcdf.netcdf_file(outfile, 'r') as f:
        z = np.array(f.variables['z'][:]).squeeze()
        T = np.array(f.variables['T'][:]).squeeze()
        act_m0 = np.array(f.variables['act_m0'][:]).squeeze()
        act_m1 = np.array(f.variables['act_m1'][:]).squeeze()
        act_m2 = np.array(f.variables['act_m2'][:]).squeeze()
        ice_m0 = np.array(f.variables['ice_m0'][:]).squeeze()
        ice_m1 = np.array(f.variables['ice_m1'][:]).squeeze()
        ice_m2 = np.array(f.variables['ice_m2'][:]).squeeze()

        ice_mix_ratio = np.array(f.variables['ice_mix_ratio'][:]).squeeze()
        liq_mix_ratio = np.array(f.variables['liq_m3'][:]).squeeze() * 4/3 * np.pi * common.rho_w
        ice_conc = np.array(f.variables['ice_m0'][:]).squeeze()  # 1/kg
        act_conc = np.array(f.variables['act_m0'][:]).squeeze()
        ice_r = np.where(ice_conc > 0, np.array(f.variables['ice_m1'][:]).squeeze() / np.array(f.variables['ice_m0'][:]).squeeze(), 0)
        act_r = np.where(act_conc > 0, np.array(f.variables['act_m1'][:]).squeeze() / np.array(f.variables['act_m0'][:]).squeeze(), 0)
        sstp_cond_mean = np.array(f.variables['sstp_cond_mean'][:]) if 'sstp_cond_mean' in f.variables else np.zeros(z.shape)
        sstp_dep_mean = np.array(f.variables['sstp_dep_mean'][:]) if 'sstp_dep_mean' in f.variables else np.zeros(z.shape)

        std_dev_liq = np.sqrt(np.where(act_m0 > 0, 
                           act_m2 / act_m0 - (act_m1 / act_m0)**2, 
                           0))
        std_dev_ice = np.sqrt(np.where(ice_m0 > 0, 
                           ice_m2 / ice_m0 - (ice_m1 / ice_m0)**2, 
                           0))
        step_cond_walltime_ms = np.array(f.variables['step_cond_walltime_ms'][:]).squeeze()

    os.remove(outfile)
    return z/1000, ice_mix_ratio*1e3, liq_mix_ratio*1e3, ice_conc/1e6, act_conc/1e6, ice_r*1e6, act_r*1e6, sstp_cond_mean, sstp_dep_mean, std_dev_liq*1e6, std_dev_ice*1e6, step_cond_walltime_ms, T



def make_figure(aerosol, w_max, z_max):

    fig, ax = plt.subplots(2, 5, figsize=(17.5, 10.0), sharey=True, squeeze=True)

    for adaptive, epsilon, sstp_act, sstp_cond in [(False, None, None, 10), (True, None, 10, 1), (True, 1e-1, None, 10),(True, 1e-2, None, 10), (True, 1e-3, None, 10)]:
    #for adaptive, epsilon, sstp_act, sstp_cond in [(True, 1e-3, None, 10)]:
    
        z, ice_mix_ratio, liq_mix_ratio, ice_conc, act_conc, ice_r, act_r, sstp_cond_mean, sstp_dep_mean, std_dev_liq, std_dev_ice, step_cond_walltime_ms, T = run_scheme(aerosol, w_max, z_max, adaptive, epsilon, sstp_act, sstp_cond)           

        walltime_mean = np.nanmean(step_cond_walltime_ms)

        if not adaptive:
            c = 'darkgrey'
            s = '-'
            l = '10 substeps, '+str(round(walltime_mean,1))+' ms'
            lw = 3
        elif sstp_act == 10:
            c = 'gold'
            s = ':'
            l = '10 activation substeps, '+str(round(walltime_mean,1))+' ms'
            lw = 2
        else:
            s = ':'
            l = 'adaptive, $\epsilon$ = '+str(epsilon)+', '+str(round(walltime_mean,1))+' ms'
            lw = 2
            if epsilon == 1e-3:
                c = 'blue'
            elif epsilon == 1e-2:
                c = 'seagreen'
            else:
                c = 'violet'

        ax[1,0].plot(ice_mix_ratio, z, color=c, label=l, linestyle=s, linewidth = lw)
        ax[0,0].plot(liq_mix_ratio, z, color=c, label=l, linestyle=s, linewidth = lw)
        ax[1,1].plot(ice_conc, z, color=c, label=l, linestyle=s, linewidth = lw)
        ax[0,1].plot(act_conc, z, color=c, label=l, linestyle=s, linewidth = lw)
        ax[1,2].plot(ice_r, z, color=c, linestyle=s, linewidth = lw)
        ax[0,2].plot(act_r, z, color=c, linestyle=s, linewidth = lw)
        ax[1,3].plot(std_dev_ice, z, color=c, label=l, linestyle=s, linewidth = lw)
        ax[0,3].plot(std_dev_liq, z, color=c, label=l, linestyle=s, linewidth = lw)
        #ax[1,3].plot(T, z, color=c, label=l, linestyle=s, linewidth = lw)
        #ax[0,3].plot(T, z, color=c, label=l, linestyle=s, linewidth = lw)
        if adaptive:
            ax[1,4].plot(sstp_dep_mean, z, color=c, linestyle=s, linewidth = lw)
            ax[0,4].plot(sstp_cond_mean, z, color=c, linestyle=s, linewidth = lw)


        for var, axis in [
                          (liq_mix_ratio, ax[0,0]), (ice_mix_ratio, ax[1,0]),
                          (act_conc, ax[0,1]), (ice_conc, ax[1,1]),
                          (act_r, ax[0,2]), (ice_r, ax[1,2]),
                          (std_dev_liq, ax[0,3]), (std_dev_ice, ax[1,3]),
                          (sstp_cond_mean, ax[0,4]), (sstp_dep_mean, ax[1,4])
                           ]:
            if not adaptive and (var is sstp_cond_mean or var is sstp_dep_mean):
                continue
            
            axis.annotate(
                '',
                xy=(var[len(var)//4+8], z[len(z)//4+8]),    # Grot strzałki
                xytext=(var[len(var)//4+6], z[len(z)//4+6]), # Początek strzałki
                arrowprops=dict(
                    arrowstyle='simple',
                    color=c,
                    linewidth=2
                )
            )
            axis.annotate(
                '',
                xy=(var[3*len(var)//5+8], z[3*len(z)//5+8]),    # Grot strzałki
                xytext=(var[3*len(var)//5+6], z[3*len(z)//5+6]), # Początek strzałki
                arrowprops=dict(
                    arrowstyle='simple',
                    color=c,
                    linewidth=2
                )
            )

    ax[0,0].set_ylabel('z [km]')
    ax[1,0].set_ylabel('z [km]')
    ax[0,0].set_xlabel('liquid mixing ratio [g/kg]')
    ax[1,0].set_xlabel('ice mixing ratio [g/kg]')
    ax[0,1].set_xlabel('liquid conc. [1/mg]')
    ax[1,1].set_xlabel('ice conc. [1/mg]')
    ax[0,2].set_xlabel("liquid mean radius [$\mu$m]")
    ax[1,2].set_xlabel("ice mean radius [$\mu$m]")
    ax[0,3].set_xlabel('liquid std. deviation [$\mu$m]')
    ax[1,3].set_xlabel('ice std. deviation [$\mu$m]')
    ax[0,4].set_xlabel("liquid sstp mean")
    ax[1,4].set_xlabel("ice sstp mean")

    ax[0,4].set_xlim(0,10.5)
    ax[1,4].set_xlim(0,10.5)

    handles, labels = ax[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="center right", bbox_to_anchor=(0.22, 0.5))
    fig.tight_layout(rect=[0.21, 0, 1, 0.95])
    aerosol_str = "pristine" if aerosol==pristine else "polluted"
    out_png = "plots/outputs/adaptive/test_adaptive_WBF_w_"+str(w_max)+"_"+aerosol_str+".pdf"
    plt.suptitle('w = '+str(w_max)+' m/s, '+aerosol_str)
    plt.savefig(out_png, dpi=200)

    return fig

for w_max,z_max in zip(w_list, z_max_list):
    for aerosol in aerosol_list:
        make_figure(aerosol, w_max, z_max)
        plt.show()