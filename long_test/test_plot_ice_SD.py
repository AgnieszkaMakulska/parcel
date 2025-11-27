import sys, os, subprocess
sys.path.insert(0, "../../")
sys.path.insert(0, "./")
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from parcel import parcel
from scipy.io import netcdf
from libcloudphxx import common

def plot_profiles(fnc, output_folder="plots/outputs/"):

    plt.figure(1, figsize=(15,5))
    plots    = []
    for i in range(2):
        plots.append(plt.subplot(1,2,i+1))
    plots[0].set_xlabel('mixing ratio [g/kg]')
    plots[1].set_xlabel('T [K]')
    for ax in plots:
        ax.set_ylabel('z [m]')

    z = fnc.variables["z"][:]
    r_v = fnc.variables["r_v"][:] * 1000  #g/kg
    r_liq = np.array([i[0] for i in fnc.variables['liq_m3'][:]]) *4/3 * np.pi * common.rho_w * 1000 #g/kg
    r_ice = fnc.variables["ice_mass"][:] * 1000 #g/kg

    plots[0].plot(r_v + r_liq + r_ice, z)
    plots[0].plot(r_v, z)
    plots[0].plot(r_liq, z)
    plots[0].plot(r_ice, z)
    plots[0].legend(['r_tot', 'r_v', 'r_liq', 'r_ice'], loc='best')
    plots[1].plot(fnc.variables["T"][:], z)

    if not os.path.exists(output_folder):
        subprocess.call(["mkdir", output_folder])
    plt.savefig(os.path.join(output_folder, "plot_profiles_ice_SD.png"))

def test_plot_ice_SD():   
    outfile = "onesim_plot.nc"
    parcel(dt=1.,w=1.,sd_conc=100,
           z_max = 5000.0,
           T_0 = 263.0,
           RH_0 = 1.,
           scheme = "lgrngn",
           ice_switch=True,
           ice_nucl=True,
           time_dep_ice_nucl=False,
           aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "rd_insol": 0.5e-6, "mean_r": [0.02e-6], "gstdev": [1.4], "n_tot": [60.0e6]}}', 
           outfreq = 10, 
           out_bin= '{"liq": {"rght": 1, "moms": [0,3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 5e-20}}',
           outfile=outfile)
    fnc = netcdf.netcdf_file(outfile)
    plot_profiles(fnc)
    fnc.close()
    subprocess.call(["rm", outfile])