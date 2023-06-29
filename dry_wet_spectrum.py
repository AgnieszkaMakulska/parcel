import sys, os, subprocess
sys.path.insert(0, "../../")
sys.path.insert(0, "./")
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from parcel import parcel
from scipy.io import netcdf

sstp=1
sd_conc=2000
w=2.
#dt_list=[0.1,0.4,0.7,0.8,0.9,1.]
dt_list=[1.]

def read(fnc):
    dry_moment0=[i[0] for i in fnc.variables['dry_m0'][:]]
    dry_spectrum=[i for i in fnc.variables['dry_spectrum_m0'][:]]
    rd=[i for i in fnc.variables['dry_spectrum_r_dry']]
    wet_moment0=[i[0] for i in fnc.variables['wet_m0'][:]]
    wet_spectrum=[i for i in fnc.variables['wet_spectrum_m0'][:]]
    rw=[i for i in fnc.variables['wet_spectrum_r_wet']]
    activated_moment0=[i for i in fnc.variables['rw_ge_rc_mom0'][:]]
    activated_moment3=[i for i in fnc.variables['rw_ge_rc_mom3'][:]]
    return [dry_moment0,dry_spectrum,rd, wet_moment0, wet_spectrum,rw, activated_moment0, activated_moment3]

def dry_spectrum(data, output_folder="../outputs"):
    plt.figure(1, figsize=(7,5))
    for i in range(len(data)):
        dry_moment0=data[i][0]
        dry_spectrum = np.array(data[i][1])
        r=np.array(data[i][2])
        n_tot_dry=np.multiply(dry_moment0,1e-6)
        plt.scatter(r*1e6,dry_spectrum[0]*1e-6,marker='.')
        plt.xscale('log')
        plt.xlabel('$r_d$ [$\mu$m]')
        plt.ylabel('n [1/mg]')
        plt.xlim(0.01,0.6)
    if not os.path.exists(output_folder):
        subprocess.call(["mkdir", output_folder])
    plt.title('Dry radius spectrum cumulus',fontsize=10)
    plt.savefig(os.path.join(output_folder, "dry_spectrum_.png"),dpi=300)


def wet_spectrum(data, output_folder="../outputs"):
    plt.figure(1, figsize=(7,5))
    for i in range(len(data)):
        wet_moment0=data[i][3]
        wet_spectrum = np.array(data[i][4])
        rw=np.array(data[i][5])
        n_tot_wet=np.multiply(wet_moment0,1e-6)
        plt.scatter(rw*1e6,wet_spectrum[0]*1e-6,marker='.')
        #plt.scatter(rw*1e6,wet_spectrum[int(40/(w*dt_list[i]))]*1e-6,marker='.')
        plt.scatter(rw*1e6,wet_spectrum[-1]*1e-6,marker='.') #[int(200/(w*dt_list[i]))]*1e-6,marker='.')
        plt.legend(['z=0 m','z=200 m'])
        plt.xscale('log')
        plt.xlabel('$r_w$ [$\mu$m]')
        plt.ylabel('n [1/mg]')
    if not os.path.exists(output_folder):

        subprocess.call(["mkdir", output_folder])
    plt.title('Wet radius spectrum stratocumulus, dt='+str(dt_list[i])+' s, w='+str(w)+' m/s',fontsize=10)
    plt.savefig(os.path.join(output_folder, "wet_spectrum_.png"),dpi=300)


def main():
    data=[]
    for dt in dt_list:
        outfile = "dry_spectrum.nc"
        parcel(dt=dt, w=w, sd_conc=sd_conc, 
                #aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.04e-6], "gstdev": [1.7], "n_tot": [300.0e6]}}', #default
                #aerosol = '{"ammonium_sulfate": {"kappa": 1.28, "mean_r": [0.02e-6], "gstdev": [1.4], "n_tot": [300.0e6]}}',
                #aerosol = '{"ammonium": {"kappa": 0.61, "mean_r": [0.06e-6,0.011e-6], "gstdev": [1.7,1.2], "n_tot": [65.0e6,125.0e6]} }', #stratocumulus
                aerosol = '{ "ammonium": {"kappa": 0.61, "mean_r": [0.03e-6,0.14e-6], "gstdev": [1.28,1.75], "n_tot": [90.0e6,15.0e6]} }', #cumulus
                outfreq = 1, outfile=outfile,
                #out_bin= '{"wet": {"rght": 15e-6, "moms": [0], "drwt": "wet", "nbin": 1, "lnli": "log", "left": 1e-06},"wet_spectrum": {"rght": 15e-6, "moms": [0], "drwt": "wet", "nbin": 100, "lnli": "log", "left": 1e-06},  "dry": {"rght": 0.3e-6,  "moms": [0], "drwt": "dry", "nbin": 1, "lnli": "log", "left": 1e-10},  "dry_spectrum": {"rght": 0.3e-6,  "moms": [0], "drwt": "dry", "nbin": 100, "lnli": "log", "left": 1e-10} }',
                out_bin= '{"wet_spectrum": {"rght": 13e-6, "moms": [0], "drwt": "wet", "nbin": 200, "lnli": "log", "left": 0.01e-6}, "wet": {"rght": 10e-6, "moms": [0], "drwt": "wet", "nbin": 1, "lnli": "log", "left": 0.01e-6}, "dry": {"rght": 0.5e-6,  "moms": [0], "drwt": "dry", "nbin": 1, "lnli": "log", "left": 0.005e-6},  "dry_spectrum": {"rght": 0.6e-6,  "moms": [0], "drwt": "dry", "nbin": 200, "lnli": "log", "left": 0.01e-6} }',
                sstp_cond = sstp #number of condensation substeps
                )
        fnc=(netcdf.netcdf_file(outfile))
        data.append(read(fnc))
        fnc.close()
        #subprocess.call(["rm", outfile])

    dry_spectrum(data)
    #wet_spectrum(data)
    

if __name__ == '__main__':
    main()
