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
plt.rcParams['font.size'] = 20
plt.rcParams['lines.linewidth'] = 3

sd_conc=1000
w=2.
sstp=1
dt = 1

def read(fnc):
    z = fnc.variables["z"][:]
    rhod=fnc.variables['rhod'][:]
    moment0=[i[0] for i in fnc.variables['cloud_m0'][:]]
    moment3=[i[0] for i in fnc.variables['cloud_m3'][:]]
    dry_moment0=[i[0] for i in fnc.variables['dry_m0'][:]]
    activated_moment0=[i for i in fnc.variables['rw_ge_rc_mom0'][:]]
    activated_moment3=[i for i in fnc.variables['rw_ge_rc_mom3'][:]]
    RH = fnc.variables["RH"][:]  
    return [z,rhod,RH,moment0,moment3,dry_moment0,activated_moment0,activated_moment3]

def profiles(data, output_folder="outputs"):
    fig = plt.figure(1, figsize=(15,10))
    ax = plt.subplot(111)
    ax.set_ylabel('z')    
    for i in range(len(data)):
        z = data[i][0]
        rhod=data[i][1]
        RH = data[i][2]
        moment0=data[i][3]
        moment3=data[i][4]
        dry_moment0=data[i][5]
        activated_moment0=data[i][6]
        activated_moment3=data[i][7]
        #n_tot_dry=dry_moment0*1e-6 # 1/mg
        #n_tot=moment0*1e-6 #1/mg
        #n_activated=activated_moment0*1e-6 #1/mg
        #n_tot_dry=np.multiply(dry_moment0,rhod)*1e-6 #1/cm^3
        n_tot= np.multiply(moment0,rhod)*1e-6 #1/cm-3 
        #n_activated=np.multiply(activated_moment0,rhod)*1e-6 #1/mg
        rv3=np.divide(moment3,moment0) #m^3
        rv=np.cbrt(rv3)*1e6 #\mu m
        LWC=np.multiply(moment0,rv3)*4*np.pi/3 *1e6 #g/kg
        ax.plot(n_tot/100,z,label='concentration x 100 cm$^{-3}$',c='mediumblue')
        ax.plot((RH-1)*100, z, label='supersaturation [%]',c='magenta')
        ax.plot(rv/10,z, label='mean volume radius x 10 $\mu$m',c='limegreen')
        ax.plot(LWC,z,label='liquid water content [g/kg]',c='orange')
        ax.set_xlim(-0.1,1.1)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.5, box.height])
        ax.legend(loc='center left', bbox_to_anchor=(1.1, 0.5))
    if not os.path.exists(output_folder):
        subprocess.call(["mkdir", output_folder])
    plt.suptitle('w='+str(w)+' m/s')
    plt.savefig(os.path.join(output_folder, "animation.png"),dpi=300)



def main():
    data=[]
    outfile = "concentration_saturation.nc"
    parcel(dt=dt, w=w, sd_conc=sd_conc, 
            #aerosol = '{"ammonium_sulfate": {"kappa": 1.28, "mean_r": [0.02e-6], "gstdev": [1.4], "n_tot": [300.0e6]}}', #default
            #aerosol = '{"ammonium": {"kappa": 0.61, "mean_r": [0.06e-6,0.011e-6], "gstdev": [1.7,1.2], "n_tot": [65.0e6,125.0e6]} }', #stratocumulus
            aerosol = '{ "ammonium": {"kappa": 0.61, "mean_r": [0.03e-6,0.14e-6], "gstdev": [1.28,1.75], "n_tot": [90.0e6,15.0e6]} }', #cumulus
            outfreq = 1, outfile=outfile,
            out_bin= '{"cloud": {"rght": 2.5e-05, "moms": [0,3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 5e-07},  "dry": {"rght": 1e-4,  "moms": [0], "drwt": "dry", "nbin": 1, "lnli": "lin", "left": 1e-10}}',
            sstp_cond = sstp #number of condensation substeps
            )
    fnc=(netcdf.netcdf_file(outfile))
    data.append(read(fnc))
    fnc.close()
    subprocess.call(["rm", outfile])

    profiles(data)

if __name__ == '__main__':
    main()
