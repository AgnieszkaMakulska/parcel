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


sd_conc=1000
w=2.
dt=1.
sstp_list=[1,2,5,10]
#sstp=1
#dt_list=[0.1,0.4,0.6,0.8,1.]
#dt=0.1
#w_list=[0.5,1.,2.]

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

def profiles_dt(data, output_folder="../outputs"):
    plt.figure(1, figsize=(22,10))
    plots = []
    for i in range(4):
        plots.append(plt.subplot(1,4,i+1))
        plots[i].set_ylabel('z [m]')
        plots[i].grid()
        plots[i].set_ylim(0,200)
    plots[0].set_xlabel('$n_{tot}~[1/mg]$')
    plots[0].set_title('Concentration of cloud droplets',fontsize=15)
    plots[1].set_xlabel('RH')
    plots[1].set_title('Relative humidity',fontsize=15)
    plots[2].set_xlabel(r'$r_v~[\mu m]$')
    plots[2].set_title('Mean volume droplet radius',fontsize=15)
    plots[3].set_xlabel('$q_c$ [g/kg]')
    plots[3].set_title('Specific liquid water content',fontsize=15)
    legend0=[]
    legend1=[]
    legend2=[]
    legend3=[]
    for i in range(len(data)):
        z = data[i][0]
        rhod=data[i][1]
        RH = data[i][2]
        moment0=data[i][3]
        moment3=data[i][4]
        dry_moment0=data[i][5]
        activated_moment0=data[i][6]
        activated_moment3=data[i][7]
        n_tot_dry=np.multiply(dry_moment0,1e-6) # 1/mg
        n_tot=np.multiply(moment0,1e-6) #1/mg
        n_activated=np.multiply(activated_moment0,1e-6) #1/mg
        #n_tot_dry=np.multiply(dry_moment0,rhod)*1e-6 #1/cm^3
        #n_tot= np.multiply(moment0,rhod)*1e-6 #1/cm-3 
        rv3=np.divide(activated_moment3,activated_moment0) #m^3
        rv=np.cbrt(rv3)*1e6 #\mu m
        LWC=np.multiply(activated_moment0,rv3)*4*np.pi/3 *1e6 #g/kg
        plots[0].plot(n_activated ,z)
        plots[1].plot(RH, z)
        plots[2].plot(rv,z)
        plots[3].plot(LWC,z)
        #legend0.append('$dt=$'+str(dt_list[i])+'s, $n_{max}$='+str(int(round(max(n_tot)))))
        legend0.append('$dt=$'+str(dt_list[i])+'s, $n$='+str(int(round(n_activated[-1]))))
        legend1.append('$dt=$'+str(dt_list[i])+'s, $S_{max}$='+str(round(100*(max(RH)-1),2))+'%, $z_{max}$='+str(int(round(z[np.argmax(RH)])))+'m')
        legend2.append('$dt=$'+str(dt_list[i])+'s, $r_{v~max}$='+str(round(max(rv[~np.isnan(rv)]),1)))
        legend3.append('$dt=$'+str(dt_list[i])+'s, $q_{c,max}$='+str(round(max(LWC[~np.isnan(LWC)]),3)))
        #Powiększenie
        '''
        legend2.append('$dt=$'+str(dt_list[i])+'s, $r_{v~max}$='+str(round(rv[int(75/dt_list[i])],1)))
        legend3.append('$dt=$'+str(dt_list[i])+'s, $q_{c,max}$='+str(round(LWC[int(75/dt_list[i])],3)))
        '''
    plots[0].legend(legend0, loc="upper left")
    plots[1].legend(legend1, loc="upper left")
    plots[2].legend(legend2, loc="upper left")
    plots[3].legend(legend3, loc="upper left")
    plots[1].set_xlim(1.,np.max(data[0][2])+0.002)
    plots[0].set_xlim(0.4*np.max(np.multiply(data[0][6],1e-6)),1.2*np.max(np.multiply(data[-1][6],1e-6)))
    z0=z[:len(LWC[np.isnan(LWC)])+1]
    plots[3].plot(np.zeros(z0.shape), z0, color='grey')
#Powiększenie:
    '''
    plots[2].set_xlim(0.5,3.)
    plots[3].set_xlim(-0.001,0.03)
    for i in range(4):
        plots[i].set_ylim(25,75)
    plots[2].set_xlim(0.,5.)
    plots[3].set_xlim(-0.01,0.08)
    '''
    if not os.path.exists(output_folder):
        subprocess.call(["mkdir", output_folder])
    #plt.suptitle('w='+str(w)+' m/s, sd_conc='+str(sd_conc)+', $\kappa$=1.28, $n_{aerosol}$= '+str(round(n_tot_dry[0]))+r'$~\frac{1}{\mathrm{mg}}$',fontsize=15)
    plt.suptitle('stratocumulus, w='+str(w)+' m/s, sstp='+str(sstp),fontsize=15)
    plt.savefig(os.path.join(output_folder, "concentration_saturation_.png"),dpi=300)


def profiles_sstp(data, output_folder="../outputs"):
    plt.figure(1, figsize=(22,10))
    plots = []
    for i in range(4):
        plots.append(plt.subplot(1,4,i+1))
        plots[i].set_ylabel('z [m]')
        plots[i].grid()
        plots[i].set_ylim(0,200)
    plots[0].set_xlabel('$n_{tot}~[1/mg]$')
    plots[0].set_title('Concentration of cloud droplets',fontsize=15)
    plots[1].set_xlabel('RH')
    plots[1].set_title('Relative humidity',fontsize=15)
    plots[2].set_xlabel(r'$r_v~[\mu m]$')
    plots[2].set_title('Mean volume droplet radius',fontsize=15)
    plots[3].set_xlabel('$q_c$ [g/kg]')
    plots[3].set_title('Specific liquid water content',fontsize=15)
    legend0=[]
    legend1=[]
    legend2=[]
    legend3=[]
    for i in range(len(data)):
        z = data[i][0]
        rhod=data[i][1]
        RH = data[i][2]
        moment0=data[i][3]
        moment3=data[i][4]
        dry_moment0=data[i][5]
        activated_moment0=data[i][6]
        activated_moment3=data[i][7]
        n_tot_dry=np.multiply(dry_moment0,1e-6) # 1/mg
        n_tot=np.multiply(moment0,1e-6) #1/mg
        n_activated=np.multiply(activated_moment0,1e-6) #1/mg
        #n_tot_dry=np.multiply(dry_moment0,rhod)*1e-6 #1/cm^3
        #n_tot= np.multiply(moment0,rhod)*1e-6 #1/cm-3 
        rv3=np.divide(activated_moment3,activated_moment0) #m^3
        rv=np.cbrt(rv3)*1e6 #\mu m
        LWC=np.multiply(activated_moment0,rv3)*4*np.pi/3 *1e6 #g/kg
        plots[0].plot(n_activated,z)
        plots[1].plot(RH, z)
        plots[2].plot(rv,z)
        plots[3].plot(LWC,z)
        legend0.append('$sstp=$'+str(sstp_list[i])+', $n$='+str(int(round(n_activated[-1]))))
        legend1.append('$sstp=$'+str(sstp_list[i])+', $S_{max}$='+str(round(100*(max(RH)-1),2))+'%, $z_{max}$='+str(int(round(z[np.argmax(RH)])))+'m')
        legend2.append('$sstp=$'+str(sstp_list[i])+', $r_{v~max}$='+str(round(max(rv[~np.isnan(rv)]),1)))
        legend3.append('$sstp=$'+str(sstp_list[i])+', $q_{c,max}$='+str(round(max(LWC[~np.isnan(LWC)]),3)))
        #Powiększenie
        '''
        legend2.append('$dt=$'+str(dt_list[i])+'s, $r_{v~max}$='+str(round(rv[int(75/dt_list[i])],1)))
        legend3.append('$dt=$'+str(dt_list[i])+'s, $q_{c,max}$='+str(round(LWC[int(75/dt_list[i])],3)))
        '''
    plots[0].legend(legend0, loc="upper left")
    plots[1].legend(legend1, loc="upper left")
    plots[2].legend(legend2, loc="upper left")
    plots[3].legend(legend3, loc="upper left")
    plots[1].set_xlim(1.,np.max(data[-1][2])+0.002)
    plots[0].set_xlim(0.4*np.max(np.multiply(data[-1][6],1e-6)),1.2*np.max(np.multiply(data[0][6],1e-6)))
    z0=z[:len(LWC[np.isnan(LWC)])+1]
    plots[3].plot(np.zeros(z0.shape),z0,color='grey')
#Powiększenie:
    '''
    plots[2].set_xlim(0.5,3.)
    plots[3].set_xlim(-0.001,0.03)
    for i in range(4):
        plots[i].set_ylim(25,75)
    plots[2].set_xlim(0.,5.)
    plots[3].set_xlim(-0.01,0.08)
    '''
    if not os.path.exists(output_folder):
        subprocess.call(["mkdir", output_folder])
    plt.suptitle('cumulus, w='+str(w)+' m/s, dt='+str(dt)+'s',fontsize=15)
    plt.savefig(os.path.join(output_folder, "concentration_saturation.png"),dpi=300)


def profiles_v(data, output_folder="../outputs"):

    plt.figure(1, figsize=(22,10))
    plots = []
    for i in range(4):
        plots.append(plt.subplot(1,4,i+1))
        plots[i].set_ylabel('z [m]')
        plots[i].grid()
        plots[i].set_ylim(0,200)
    #plots[0].set_xlabel('$n_{tot}~[1/cm^3]$')
    plots[0].set_xlabel('$n_{tot}~[1/mg]$')
    plots[0].set_title('Concentration of cloud droplets',fontsize=15)
    plots[1].set_xlabel('RH')
    plots[1].set_title('Relative humidity',fontsize=15)
    plots[2].set_xlabel(r'$r_v~[\mu m]$')
    plots[2].set_title('Mean volume droplet radius',fontsize=15)
    plots[3].set_xlabel('$q_c$ [g/kg]')
    plots[3].set_title('Specific liquid water content',fontsize=15)

    legend0=[]
    legend1=[]
    legend2=[]
    legend3=[]
    for i in range(len(data)):
        z = data[i][0]
        rhod=data[i][1]
        RH = data[i][2]
        moment0=data[i][3]
        moment3=data[i][4]
        dry_moment0=data[i][5]
        activated_moment0=data[i][6]
        activated_moment3=data[i][7]
        n_tot_dry=np.multiply(dry_moment0,1e-6) # 1/mg
        n_tot=np.multiply(moment0,1e-6) #1/mg
        n_activated=np.multiply(activated_moment0,1e-6) #1/mg
        #n_tot_dry=np.multiply(dry_moment0,rhod)*1e-6 #1/cm^3
        #n_tot= np.multiply(moment0,rhod)*1e-6 #1/cm-3 
        rv3=np.divide(activated_moment3,activated_moment0) #m^3
        rv=np.cbrt(rv3)*1e6 #\mu m
        LWC=np.multiply(activated_moment0,rv3)*4*np.pi/3 *1e6 #g/kg
        #plots[0].plot(n_tot ,z)
        plots[0].plot(n_activated ,z)
        plots[1].plot(RH, z)
        plots[2].plot(rv,z)
        plots[3].plot(LWC,z)
        legend0.append('w='+str(w_list[i])+'m/s, $n$='+str(int(round(n_activated[-1]))))
        legend1.append('w='+str(w_list[i])+'m/s, $S_{max}$='+str(round(100*(max(RH)-1),2))+'%, $z_{max}$='+str(int(round(z[np.argmax(RH)])))+'m')
        legend2.append('w='+str(w_list[i])+'m/s, $r_{v~max}$='+str(round(max(rv[~np.isnan(rv)]),1)))
        legend3.append('w='+str(w_list[i])+'m/s, $q_{c,max}$='+str(round(max(LWC[~np.isnan(LWC)]),3)))
    plots[0].legend(legend0, loc="upper left")
    plots[1].legend(legend1, loc="upper left")
    plots[2].legend(legend2, loc="upper left")
    plots[3].legend(legend3, loc="upper left")
    #plots[1].set_xlim(max(RH)-0.002,max(RH)+0.002)
    #plots[0].set_xlim(max(n_tot)*0.7,max(n_tot)*1.05)
    z0=z[:len(LWC[np.isnan(LWC)])+1]
    plots[3].plot(np.zeros(z0.shape), z0, color='grey')
   
    if not os.path.exists(output_folder):
        subprocess.call(["mkdir", output_folder])
    plt.suptitle('dt='+str(dt)+' s, sd_conc='+str(sd_conc),fontsize=15)
    plt.savefig(os.path.join(output_folder, "concentration_saturation_.png"),dpi=300)


def main():
    data=[]
    #for dt in dt_list:
    #for w in w_list:
    for sstp in sstp_list:
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

    #profiles_dt(data)
    #profiles_v(data)
    profiles_sstp(data)

    

if __name__ == '__main__':
    main()
