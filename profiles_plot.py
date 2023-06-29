import sys, os, subprocess
sys.path.insert(0, "../../")
sys.path.insert(0, "./")

from scipy.io import netcdf
import numpy as np

from parcel import parcel

sd_conc=100
w=1.
dt=0.1

def plot_profiles(fnc, output_folder="../outputs"):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties

    # ... plotting the results ...
    plt.figure(1, figsize=(18,10))
    plots    = []
    legend_l = []
    for i in range(6):
        plots.append(plt.subplot(2,3,i+1))
    

    plots[0].set_xlabel('p [hPa]')

    plots[1].ticklabel_format(useOffset=False) 
    plots[1].set_xlabel('th_d [K]')
    plots[2].set_xlabel('T [K]')
    plots[3].set_xlabel('rho_d [kg/m3]')  
    plots[4].set_xlabel('rv [g/kg]')
    plots[5].set_xlabel('RH')

    for ax in plots:
        ax.set_ylabel('z [m]')

    z = fnc.variables["z"][:]
    plots[0].plot(fnc.variables["p"][:] / 100.   , z)
    plots[1].plot(fnc.variables["th_d"][:]       , z)
    plots[2].plot(fnc.variables["T"][:]          , z)
    plots[3].plot(fnc.variables["rhod"][:]      , z)
    plots[4].plot(fnc.variables["r_v"][:] * 1000 , z)
    plots[5].plot(
        fnc.variables["RH"][:]                     , z, 
        [fnc.variables["RH"][:].max()] * z.shape[0], z
        )
    
   
    if not os.path.exists(output_folder):
        subprocess.call(["mkdir", output_folder])
    plt.suptitle('dt='+str(dt)+' s, w='+str(w)+' m/s, sd_conc='+str(sd_conc)+', $\kappa$=1.28',fontsize=18)
    plt.savefig(os.path.join(output_folder, "plot_profiles_.png"))

def main():
    # running parcel model for different ways to solve for pressure  ...            
    outfile = "onesim_plot.nc"
    parcel(dt=0.1,w=1.,sd_conc=100,aerosol = '{"ammonium_sulfate": {"kappa": 1.28, "mean_r": [0.02e-6], "gstdev": [1.4], "n_tot": [60.0e6]}}', outfreq = 10, outfile=outfile)
    fnc = netcdf.netcdf_file(outfile)
    plot_profiles(fnc)
    fnc.close()
    subprocess.call(["rm", outfile])
    

if __name__ == '__main__':
    main()
