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
plt.rcParams['font.size'] = 12
import csv

sd_conc=1000
sstp=1
dt=0.1
w_list=[1.,2.]
#dt=1.
#w_list=[1.]

def read(fnc):
    z = fnc.variables["z"][:]
    rhod=fnc.variables['rhod'][:]
    moment0=[i[0] for i in fnc.variables['cloud_m0'][:]]
    moment3=[i[0] for i in fnc.variables['cloud_m3'][:]]
    dry_moment0=[i[0] for i in fnc.variables['dry_m0'][:]]
    activated_moment0=[i for i in fnc.variables['rw_ge_rc_mom0'][:]]
    activated_moment3=[i for i in fnc.variables['rw_ge_rc_mom3'][:]]
    RH = fnc.variables["RH"][:]  
    T= fnc.variables["T"][:]  
    p= fnc.variables["p"][:]
    return [z,rhod,RH,T,p,activated_moment0,activated_moment3, moment0, moment3]

def data_v(data, output_folder="../outputs"):
    
    for i in range(len(data)):
        z = data[i][0]
        rhod=data[i][1]
        RH = data[i][2]
        T= data[i][3]
        p= data[i][4]
        activated_moment0=data[i][5]
        activated_moment3=data[i][6]
        moment0=data[i][7]
        moment3=data[i][8]
        #n_tot_dry=np.multiply(dry_moment0,1e-6) # 1/mg
        #n_tot=np.multiply(moment0,1e-6) #1/mg
        #n_activated=np.multiply(activated_moment0,1e-6) #1/mg
        #n_tot_dry=np.multiply(dry_moment0,rhod)*1e-6 #1/cm^3

        n_tot= np.multiply(activated_moment0,rhod)*1e-6 #1/cm-3 
        rv3=np.divide(activated_moment3,activated_moment0) #m^3
        rv=np.cbrt(rv3)*1e6 #\mu m
        LWC=np.multiply(moment0,np.divide(moment3,moment0))*4*np.pi/3 *1e6 #g/kg

        #n_tot= np.multiply(moment0,rhod)*1e-6 #1/cm-3 
        #rv3=np.divide(moment3,moment0) #m^3
        #rv=np.cbrt(rv3)*1e6 #\mu m
        #LWC=np.multiply(moment0,rv3)*4*np.pi/3 *1e6 #g/kg

        header =['z [m]','N [cm^-3]','RH','r_v [um]','q_c [g/kg]','T [K]','p [Pa]']
        rows=np.transpose(np.array([z,n_tot,RH,rv,LWC,T,p]))
   
        if not os.path.exists(output_folder):
            subprocess.call(["mkdir", output_folder])
        name='w='+str(int(w_list[i]))
        with open(os.path.join(output_folder, name), 'w') as f:  
            writer = csv.writer(f)
            writer.writerow(header)
            writer.writerows(rows)


def main():
    data=[]
    for w in w_list:
        outfile = "data.nc"
        parcel(dt=dt, w=w, sd_conc=sd_conc, 
                aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.04e-6], "gstdev": [1.7], "n_tot": [300.0e6]}}', #default
                #aerosol = '{"ammonium": {"kappa": 0.61, "mean_r": [0.06e-6,0.011e-6], "gstdev": [1.7,1.2], "n_tot": [65.0e6,125.0e6]} }', #stratocumulus
                #aerosol = '{ "ammonium": {"kappa": 0.61, "mean_r": [0.03e-6,0.14e-6], "gstdev": [1.28,1.75], "n_tot": [90.0e6,15.0e6]} }', #cumulus
                outfreq = 1, outfile=outfile,
                out_bin= '{"cloud": {"rght": 25e-06, "moms": [0,3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0.5e-06},  "dry": {"rght": 1e-4,  "moms": [0], "drwt": "dry", "nbin": 1, "lnli": "lin", "left": 1e-10}}',
                sstp_cond = sstp #number of condensation substeps
                )
        fnc=(netcdf.netcdf_file(outfile))
        data.append(read(fnc))
        fnc.close()
        subprocess.call(["rm", outfile])

    data_v(data)
    

if __name__ == '__main__':
    main()
