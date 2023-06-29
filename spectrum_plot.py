# This Python file uses the following encoding: utf-8
import sys
sys.path.insert(0, "../../")
sys.path.insert(0, "../")
sys.path.insert(0, "./")

from scipy.io import netcdf
import numpy as np
import pytest
import subprocess

from parcel import parcel
sd_conc=1000
w=1.
dt=0.1

def spectrum(fnc, output_folder="../outputs"):

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.font_manager import FontProperties
    import os

    ymax = 1e4
    ymin = 1e-1

    rw = fnc.variables["wradii_r_wet"][:] * 1e6
    rd = fnc.variables["dradii_r_dry"][:] * 1e6

    drw = np.empty(rw.shape) 
    drw[0] = rw[0] - 0
    drw[1:] = (rw[1:] - rw[0:-1]) * 1e6

    drd = np.empty(rd.shape)
    drd[0] = rd[0] - 0
    drd[1:] = (rd[1:] - rd[0:-1]) * 1e6

    for t in range(fnc.variables['t'].shape[0]):

        plt.clf()

        out=output_folder + 'plot_spec_' + str("%03d" % t) + '.png'
        # ... plotting the results ...
        plt.figure(1, figsize=(18,10))

        plt.ylabel('n')
        plt.xlabel('r [$\mu$m]')
        plt.suptitle('Spectrum',fontsize=15)

        plt.xscale('log')
        plt.yscale('log')
        plt.ylim(ymin,ymax)

        
        nw = fnc.variables['wradii_m0'][t,:] / drw  #n*rho/drw
        nd = fnc.variables['dradii_m0'][t,:] / drd

        plt.scatter(rw, nw)
        plt.scatter(rd, nd)
        plt.legend(['wet','dry'])
        
        #plot_rw = Gnuplot.PlotItems.Data(rw, nw, with_="fsteps", title="wet radius")
        #plot_rd = Gnuplot.PlotItems.Data(rd, nd, with_="fsteps", title="dry radius")
    
        if not os.path.exists(output_folder):
            subprocess.call(["mkdir", output_folder])
        plt.title('dt='+str(dt)+' s, w='+str(w)+' m/s, sd_conc='+str(sd_conc)+', $\kappa$=1.28',fontsize=15)
        plt.savefig(os.path.join(out))


def main():

    outfile = "test_spectrum.nc"
    out_bin = '{"wradii": {"rght": 1e-4, "left": 1e-9, "drwt": "wet", "lnli": "log", "nbin": 26, "moms": [0]},"dradii": {"rght": 1e-6, "left": 1e-9, "drwt": "dry", "lnli": "log", "nbin": 26, "moms": [0]}}'

    # run parcel run!
    parcel(dt = .5, sd_conc = 1024, outfreq = 20,  outfile = outfile, out_bin = out_bin)

    fnc = netcdf.netcdf_file(outfile, "r")

    # plotting 
    spectrum(fnc, output_folder="../outputs/")

    # cleanup
    subprocess.call(["rm", outfile])


if __name__ == '__main__':
    main()

