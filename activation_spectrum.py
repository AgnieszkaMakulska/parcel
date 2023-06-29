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
from scipy.special import erf

sstp=1
sd_conc=1000
w=1.
dt_list=[0.1,0.5,1.]
#dt_list=[1.]

def activation_formula(N,mean_r,gstdev,kappa,T,S):
    #sigma= 75.64*10**(-3)
    sigma = 0.07275 * (1-0.002 * (T-291))
    rho_l= 997
    R_v= 461.5 
    A=2*sigma/(rho_l*R_v*T)
    r_crit= ( 4*A**3 / (27*kappa*(S-1)**2) )**(1/3)
    return N/2 * (1 + erf(np.log(mean_r/r_crit) / (np.sqrt(2)*np.log(gstdev)) ) )


def read(fnc):
    z = fnc.variables["z"][:]
    rhod=fnc.variables['rhod'][:]
    RH = fnc.variables["RH"][:]  
    T = fnc.variables["T"][:] 
    moment0=[i[0] for i in fnc.variables['cloud_m0'][:]]
    dry_moment0=[i[0] for i in fnc.variables['dry_m0'][:]]
    activated_moment0=[i for i in fnc.variables['rw_ge_rc_mom0'][:]]
    critical_moment0=[i for i in fnc.variables['RH_ge_Sc_mom0'][:]]
    return [z,rhod,RH,moment0,dry_moment0,T,activated_moment0,critical_moment0]

def profiles_dt(data, output_folder="../outputs"):

    plt.figure(1, figsize=(7,5))

    legend=[]
    for i in range(len(data)):
        z = data[i][0]
        rhod=data[i][1]
        RH = data[i][2]
        moment0=data[i][3]
        dry_moment0=data[i][4]
        T = data[i][5]
        activated_moment0=data[i][6]
        critical_moment0=data[i][7]
        n_tot_dry=np.multiply(dry_moment0,1e-6) # 1/mg
        n_tot=np.multiply(moment0,1e-6) #1/mg
        n_activated=np.multiply(activated_moment0,1e-6) #1/mg
        n_crit=np.multiply(critical_moment0,1e-6) #1/mg
        #n_tot_dry=np.multiply(dry_moment0,rhod)*1e-6 #1/cm^3
        #n_tot= np.multiply(moment0,rhod)*1e-6 #1/cm-3 
        plt.plot(RH[0:np.argmax(RH)],n_crit[0:np.argmax(RH)],linestyle='--',marker='.')
        legend.append('$dt=$'+str(dt_list[i])+'s')
        #plt.plot(RH[0:np.argmax(RH)],n_activated[0:np.argmax(RH)],linestyle='--',marker='.')
        #legend.append('$dt=$'+str(dt_list[i])+'s, formula A')
        plt.xlim(1.,1.0048)
        plt.xlabel('RH')
        plt.ylabel('$n_{tot}~[1/mg]$')
    RH = data[0][2]
    T=data[0][5]
    plt.plot(RH, activation_formula(245,0.02*1e-6,1.4,1.28,T,RH)) #1/mg
    #plt.plot(RH, activation_formula(278,0.02*1e-6,1.4,1.28,T,RH)) #1/cm^3
    legend.append('theory')
    plt.legend(legend, loc="upper left")
    plt.ylim(0.,205.)

    if not os.path.exists(output_folder):
        subprocess.call(["mkdir", output_folder])
    plt.suptitle('Activation spectrum')
    plt.title('w='+str(w)+' m/s, sd_conc='+str(sd_conc)+', $\kappa$=1.28, $n_{aerosol}$= '+str(round(n_tot_dry[0]))+r'$~\frac{1}{\mathrm{mg}}$')
    plt.savefig(os.path.join(output_folder, "activation_spectrum_.png"),dpi=300)


def main():

    data=[]
    for dt in dt_list:
        outfile = "activation_spectrum.nc"
        parcel(dt=dt, w=w, sd_conc=sd_conc, 
                aerosol = '{"ammonium_sulfate": {"kappa": 1.28, "mean_r": [0.02e-6], "gstdev": [1.4], "n_tot": [300.0e6]}}',
                #aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.06e-6], "gstdev": [1.7], "n_tot": [65.0e6]} , "ammonium": {"kappa": 0.61, "mean_r": [0.011e-6], "gstdev": [1.2], "n_tot": [125.0e6]} }', #stratocumulus
                #aerosol = '{"ammonium_sulfate": {"kappa": 0.61, "mean_r": [0.03e-6], "gstdev": [1.28], "n_tot": [90.0e6]} , "ammonium": {"kappa": 0.61, "mean_r": [0.14e-6], "gstdev": [1.75], "n_tot": [15.0e6]} }', #cumulus
                outfreq = 1, outfile=outfile,
                out_bin= '{"cloud": {"rght": 2.5e-05, "moms": [0,3], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 5e-07},  "dry": {"rght": 1e-4,  "moms": [0], "drwt": "dry", "nbin": 1, "lnli": "lin", "left": 1e-10}}',
                sstp_cond = sstp #number of condensation substeps
                )
        fnc=(netcdf.netcdf_file(outfile))
        data.append(read(fnc))
        fnc.close()
        subprocess.call(["rm", outfile])

    profiles_dt(data)


    

if __name__ == '__main__':
    main()
