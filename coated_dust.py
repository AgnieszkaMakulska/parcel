import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")

import numpy as np
from parcel import parcel
from scipy.io import netcdf
from libcloudphxx import common


def mixed_aerosol(aerosol_str, epsilon, rd):

    if epsilon == 1.0:
        epsilon = 0.9999999
    if aerosol_str == "pristine":
        return f'{{"pristine":{{"kappa": 0.61, "sol_frac": 1.0, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125.0e6, 65.0e6]}}, \
            "mixed": {{"kappa": 0.61, "sol_frac": {epsilon}, "mean_r": [{rd}], "gstdev": [1.2], "n_tot": [5.0e6]}} }}'
    elif aerosol_str == "polluted":
        return f'{{"polluted":{{"kappa": 0.61, "sol_frac": 1.0, "mean_r": [0.029e-6, 0.071e-6], "gstdev": [1.36, 1.57], "n_tot": [160.0e6, 380.0e6]}}, \
            "mixed": {{"kappa": 0.61, "sol_frac": {epsilon}, "mean_r": [{rd}], "gstdev": [1.2], "n_tot": [5.0e6]}} }}'
    else:
        raise ValueError('unknown aerosol spec')

def soluble_aerosol(aerosol_str):
    if aerosol_str == "pristine":
        return '{"pristine": {"kappa": 0.61, "sol_frac": 1.0, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125.0e6, 65.0e6]}}'
    elif aerosol_str == "polluted":
        return '{"polluted": {"kappa": 0.61, "sol_frac": 1.0, "mean_r": [0.029e-6, 0.071e-6], "gstdev": [1.36, 1.57], "n_tot": [160.0e6, 380.0e6]}}'
    else:
        raise ValueError('unknown aerosol spec')



def run_scheme(aerosol, outfile, outfreq, spec=False):

    if spec == False:
        out_bin = '{"liq": {"rght": 1, "moms": [0,1,2,3,4], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 1e-20},' \
        '"aerosol": {"rght": 1, "moms": [0], "drwt": "dry", "nbin": 1, "lnli": "lin", "left": 1e-20},' \
        '"cloud": {"rght": 1, "moms": [0,1,2,3,4], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 1e-6}}'
    else:
        out_bin = '{"initial_spec": {"rght": 8e-6, "moms": [0], "drwt": "wet", "nbin": 100, "lnli": "log", "left": 0.01e-6},' \
            '"spec": {"rght": 30e-6, "moms": [0], "drwt": "wet", "nbin": 1000, "lnli": "log", "left": 1e-6}}'

    args = dict(
        p_0=90000,
        RH_0=0.97,
        T_0=283,
        aerosol = aerosol,
        w = 1,
        sd_conc = None,
        sd_const_multi=1000000,
        n_sd_max=1e7,
        dt = 1,
        z_max = 200,
        outfile = outfile,
        outfreq = outfreq,
        scheme = "lgrngn",
        out_bin = out_bin,
        sstp_cond = 10,
        ice_switch = False,
        ice_nucl = False,
        depo = False,
        backend = "gpu",
        aerosol_independent_of_rhod = True
        #large_tail = True
    )
    parcel(**args)


def read_profiles(outfile):
    with netcdf.netcdf_file(outfile, 'r') as f:
        z = np.array(f.variables['z'][:]).squeeze()
        cloud_m0 = np.array(f.variables['act_m0'][:]).squeeze()
        cloud_m1 = np.array(f.variables['act_m1'][:]).squeeze()
        cloud_m2 = np.array(f.variables['act_m2'][:]).squeeze()
        cloud_m3 = np.array(f.variables['act_m3'][:]).squeeze()
        cloud_m4 = np.array(f.variables['act_m4'][:]).squeeze()
        liq_mix_ratio = cloud_m3 * 4/3 * np.pi * common.rho_w
        conc = cloud_m0 
        mean_r = np.where(cloud_m0 > 0, cloud_m1 / cloud_m0, 0)
        variance_r = np.where(cloud_m0 > 0,
                        cloud_m2 / cloud_m0 - (cloud_m1 / cloud_m0)**2, 
                        0)
        std_dev_r = np.sqrt(np.maximum(variance_r, 0.0))
        print(np.array(f.variables['aerosol_m0'][:]).squeeze())
    return z, liq_mix_ratio*1e3, conc/1e6, mean_r*1e6, std_dev_r*1e6

def read_distr(outfile, z_distr):
    with netcdf.netcdf_file(outfile, 'r') as f:
        z = np.array(f.variables['z'][:]).squeeze()
        distr = np.array(f.variables['spec_m0'][:]).squeeze()
        radii = np.array(f.variables['spec_r_wet'][:]).squeeze()
        bin_widths = np.array(f.variables['spec_dr_wet'][:]).squeeze()
        init_distr = np.array(f.variables['initial_spec_m0'][:]).squeeze()
        init_radii = np.array(f.variables['initial_spec_r_wet'][:]).squeeze()
        init_bin_widths = np.array(f.variables['initial_spec_dr_wet'][:]).squeeze()
        initial_distr = init_distr[np.argmin(np.abs(z))]
        distr1 = distr[np.argmin(np.abs(z - z_distr))]
    return distr1/1e6, radii*1e6, bin_widths*1e6, initial_distr/1e6, init_radii*1e6, init_bin_widths*1e6