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
sd_conc=1000
w_list=[0.5,1.,2.]
#sstp=1
#dt_list=[0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.]
dt=1.
sstp_list=[1,2,3,4,5,6,7,8,9,10]

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

def n_dt(data, output_folder="../outputs"):
    plt.figure(1, figsize=(7,5))
    legend=[]
    for w_id in range(len(w_list)):
        n=[]
        for dt_id in range(len(dt_list)):
            i=w_id*len(dt_list)+dt_id
            dry_moment0=data[i][5]
            n_tot_dry=np.multiply(dry_moment0,1e-6) # 1/mg
            activated_moment0=data[i][6]
            n_activated=np.multiply(activated_moment0,1e-6) #1/mg
            n.append(n_activated[-1])
        plt.scatter(dt_list, n)
        legend.append('w='+str(w_list[w_id])+' m/s')
    plt.legend(legend)
    plt.xlabel('dt [s]')
    plt.ylabel('$n_{tot}$ [1/mg]')
    plt.grid()
    if not os.path.exists(output_folder):
        subprocess.call(["mkdir", output_folder])
    #plt.title('sd_conc='+str(sd_conc)+', $\kappa$=1.28, $n_{aerosol}$= '+str(round(n_tot_dry[0]))+r'$~\frac{1}{\mathrm{mg}}$',fontsize=15)
    plt.title('stratocumulus, sstp='+str(sstp)+', sd_conc='+str(sd_conc),fontsize=15)
    plt.savefig(os.path.join(output_folder, "n_dt.png"),dpi=300)

def n_sstp(data, output_folder="../outputs"):
    plt.figure(1, figsize=(7,5))
    legend=[]
    for w_id in range(len(w_list)):
        n=[]
        for dt_id in range(len(sstp_list)):
            i=w_id*len(sstp_list)+dt_id
            dry_moment0=data[i][5]
            n_tot_dry=np.multiply(dry_moment0,1e-6) # 1/mg
            activated_moment0=data[i][6]
            n_activated=np.multiply(activated_moment0,1e-6) #1/mg
            n.append(n_activated[-1])
        plt.scatter(sstp_list, n)
        legend.append('w='+str(w_list[w_id])+' m/s')
    plt.legend(legend)
    plt.xlabel('substeps')
    plt.ylabel('$n_{tot}$ [1/mg]')
    plt.grid()
    if not os.path.exists(output_folder):
        subprocess.call(["mkdir", output_folder])
    #plt.title('sd_conc='+str(sd_conc)+', $\kappa$=1.28, $n_{aerosol}$= '+str(round(n_tot_dry[0]))+r'$~\frac{1}{\mathrm{mg}}$',fontsize=15)
    plt.title('cumulus, dt='+str(dt)+'s, sd_conc='+str(sd_conc),fontsize=15)
    plt.savefig(os.path.join(output_folder, "n_dt.png"),dpi=300)

def main():

    data=[]
    for w in w_list:
        #for dt in dt_list:
        for sstp in sstp_list:
            outfile = "n_dt.nc"
            parcel(dt=dt, w=w, sd_conc=sd_conc, 
                    #aerosol = '{"ammonium_sulfate": {"kappa": 1.28, "mean_r": [0.02e-6], "gstdev": [1.4], "n_tot": [300.0e6]}}',
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

    #n_dt(data)
    n_sstp(data)
    
  

if __name__ == '__main__':
    main()
