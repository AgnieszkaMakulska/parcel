import numpy as np
import matplotlib.pyplot as plt
import sys, os, subprocess
sys.path.insert(0, "../../")
sys.path.insert(0, "./")
output_folder="../outputs"
#from scipy.stats import lognorm


def lognormal(r,mean_r,gstdev,n_tot):
   A=n_tot/(np.sqrt(2*np.pi)*np.log(gstdev))
   B= -(np.log(r)-np.log(mean_r))**2 / (2*np.log(gstdev)**2)
   return A/r *np.exp(B)

mean_r=0.02 #\mu m
gstdev=1.4
n_tot=245 #1/mg

r=np.linspace(0.001,0.1,100) #\mu m
n=lognormal(r,mean_r=mean_r, gstdev=gstdev, n_tot=n_tot) 
plt.plot(r,n)

#r=np.linspace(1e-9,1e-7,100) #\mu m
#n=lognormal(r,mean_r=0.02e-6, gstdev=1.6, n_tot=300) #[per cm^3]
#plt.plot(r*1e6,n/1e6)
#plt.xscale('log')
plt.xlabel('r [$\mu$m]')
plt.ylabel(r'$\frac{dN}{dr}$ $\left[ \frac{\mathrm{mg}^{-1}}{\mu \mathrm{m}} \right]$')
#plt.title('Lognormal distribution of aerosol \n mean_r='+str(mean_r)+'$\mu$m, gstdev='+str(gstdev)+', n_tot='+str(n_tot)+'cm$^{-3}$',fontsize=10)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, "aerosol.png"),dpi=300)

print('integral:')
print(np.trapz(n,x=r))
plt.show()
