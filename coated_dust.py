import sys
sys.path.insert(0, "../")
sys.path.insert(0, "./")

import numpy as np
from parcel import parcel
from scipy.io import netcdf
from libcloudphxx import common


def aerosol_spec(aerosol_str, epsilon, rd):

    if epsilon == 1.0:
        epsilon = 0.999
    if aerosol_str == "pristine":
        aerosol = f'{{"pristine":{{"kappa": 0.61, "sol_frac": 1.0, "mean_r": [0.011e-6, 0.06e-6], "gstdev": [1.2, 1.7], "n_tot": [125.0e6, 65.0e6]}}, \
            "mixed": {{"kappa": 0.61, "sol_frac": {epsilon}, "mean_r": [{rd}], "gstdev": [1.4], "n_tot": [1.0e6]}} }}'
    elif aerosol_str == "polluted":
        aerosol = f'{{"polluted":{{"kappa": 0.61, "sol_frac": 1.0, "mean_r": [0.029e-6, 0.071e-6], "gstdev": [1.36, 1.57], "n_tot": [160.0e6, 380.0e6]}}, \
            "mixed": {{"kappa": 0.61, "sol_frac": {epsilon}, "mean_r": [{rd}], "gstdev": [1.4], "n_tot": [1.0e6]}} }}'
    return aerosol


def run_scheme(aerosol, outfile, outfreq, spec=False):

    if spec == False:
        out_bin = '{"liq": {"rght": 1, "moms": [0,1,2,3,4], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0.5e-20}}'
    else:
        out_bin = '{"liq": {"rght": 1, "moms": [0,1,2,3,4], "drwt": "wet", "nbin": 1, "lnli": "lin", "left": 0.5e-20},' \
            '"initial_spec": {"rght": 3e-6, "moms": [0], "drwt": "wet", "nbin": 100, "lnli": "log", "left": 0.01e-6},' \
            '"spec": {"rght": 30e-6, "moms": [0], "drwt": "wet", "nbin": 1000, "lnli": "log", "left": 1e-6}}'

    args = dict(
        p_0=90000,
        RH_0=0.97,
        T_0=283,
        aerosol = aerosol,
        w = 1,
        sd_conc = 1000,
        #sd_const_multi=1000000,
        #n_sd_max=1e7,
        dt = 1,
        z_max = 400,
        outfile = outfile,
        outfreq = outfreq,
        scheme = "lgrngn",
        out_bin = out_bin,
        sstp_cond = 10,
        ice_switch = False,
        ice_nucl = False,
        depo = False,
        backend = "gpu"
    )
    parcel(**args)


def read_profiles(outfile):
    with netcdf.netcdf_file(outfile, 'r') as f:
        z = np.array(f.variables['z'][:]).squeeze()
        act_m0 = np.array(f.variables['act_m0'][:]).squeeze()
        act_m1 = np.array(f.variables['act_m1'][:]).squeeze()
        act_m2 = np.array(f.variables['act_m2'][:]).squeeze()
        act_m4 = np.array(f.variables['act_m4'][:]).squeeze()
        liq_m0 = np.array(f.variables['liq_m0'][:]).squeeze()
        liq_m1 = np.array(f.variables['liq_m1'][:]).squeeze()
        liq_m2 = np.array(f.variables['liq_m2'][:]).squeeze()
        liq_m3 = np.array(f.variables['liq_m3'][:]).squeeze()
        liq_m4 = np.array(f.variables['liq_m4'][:]).squeeze()
        liq_mix_ratio = liq_m3 * 4/3 * np.pi * common.rho_w
        conc = act_m0
        mean_r = np.where(act_m0 > 0, act_m1 / act_m0, 0)
        std_dev_r = np.sqrt(np.where(act_m0 > 0, 
                        act_m2 / act_m0 - (act_m1 / act_m0)**2, 
                        0))
    return z, liq_mix_ratio*1e3, conc/1e6, mean_r*1e6, std_dev_r*1e6

def read_distr(outfile):
    with netcdf.netcdf_file(outfile, 'r') as f:
        z = np.array(f.variables['z'][:]).squeeze()
        distr = np.array(f.variables['spec_m0'][:]).squeeze()
        radii = np.array(f.variables['spec_r_wet'][:]).squeeze()
        bin_widths = np.array(f.variables['spec_dr_wet'][:]).squeeze()
        init_distr = np.array(f.variables['initial_spec_m0'][:]).squeeze()
        init_radii = np.array(f.variables['initial_spec_r_wet'][:]).squeeze()
        init_bin_widths = np.array(f.variables['initial_spec_dr_wet'][:]).squeeze()
        initial_distr = init_distr[np.argmin(np.abs(z))]
        distr1 = distr[np.argmin(np.abs(z - 400))]
    return distr1/1e6, radii*1e6, bin_widths*1e6, initial_distr/1e6, init_radii*1e6, init_bin_widths*1e6