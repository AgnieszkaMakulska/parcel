#!/usr/bin/env python
import numpy as np
from libcloudphxx import lgrngn
from parcel_common import lognormal, sum_of_lognormals, _Chem_g_id, _Chem_a_id, _stats


def _micro_init(aerosol, opts, state):
  """Initialize the lagrangian microphysics scheme"""

  # lagrangian scheme options
  opts_init = lgrngn.opts_init_t()
  for opt in [
    "dt",
    "chem_rho", 
    "sstp_cond",
    "ice_switch",
    "time_dep_ice_nucl",
    "aerosol_independent_of_rhod",
    "const_p"
  ]:
    if opt in opts and opts[opt] is not None:
      setattr(opts_init, opt, opts[opt])

  if opts["rng_seed"] is not None:
      opts_init.rng_seed = int(opts["rng_seed"])

  if opts["sd_conc"] is not None and opts["sd_const_multi"] is not None:
    raise ValueError("sd_conc and sd_const_multi can't both be defined")  
  elif opts["sd_conc"] is not None:
    opts_init.sd_conc = int(opts["sd_conc"])
  elif opts["sd_const_multi"] is not None:
    opts_init.sd_const_multi = int(opts["sd_const_multi"])

  if opts["n_sd_max"] is not None:
    opts_init.n_sd_max = int(opts["n_sd_max"])

  opts_init.th_dry = True
  opts_init.const_p = False

  # read in the initial aerosol size distribution
  dry_distros = {}
  for name, dct in aerosol.items(): # loop over kappas
    lognormals = []
    for i in range(len(dct["mean_r"])):
      lognormals.append(lognormal(dct["mean_r"][i], dct["gstdev"][i], dct["n_tot"][i]))
    dry_distros[(float(dct["kappa"]), float(dct["sol_frac"]))] = sum_of_lognormals(lognormals)
  opts_init.dry_distros = dry_distros

  # better resolution for the SD tail
  if opts["large_tail"]:
      opts_init.sd_conc_large_tail = 1
      opts_init.n_sd_max = int(1e6)  # some more space for the tail SDs

  # switch off sedimentation and collisions
  opts_init.sedi_switch = False
  opts_init.coal_switch = False

  # switching on chemistry if either dissolving, dissociation or reactions are chosen
  opts_init.chem_switch = False
  if opts["chem_dsl"] or opts["chem_dsc"] or opts["chem_rct"]:
    opts_init.chem_switch = True
    opts_init.sstp_chem = opts["sstp_chem"]

  # initialisation
  backend_str = opts.get("backend", "serial")
  if backend_str is None:
    backend_str = "serial"
  backend_str = str(backend_str).lower()

  backend_map = {
    "serial": lgrngn.backend_t.serial,
    "openmp": lgrngn.backend_t.OpenMP,
    "omp": lgrngn.backend_t.OpenMP,
    "cuda": lgrngn.backend_t.CUDA,
    "gpu": lgrngn.backend_t.CUDA,
  }
  if backend_str not in backend_map:
    raise ValueError(f"Unknown lgrngn backend: {backend_str!r} (expected one of: {', '.join(sorted(backend_map))})")

  micro = lgrngn.factory(backend_map[backend_str], opts_init)

  ambient_chem = {}
  if micro.opts_init.chem_switch:
    ambient_chem = dict((v, state[k]) for k,v in _Chem_g_id.items())
  micro.init(state["th_d"], state["r_v"], state["rhod"], ambient_chem=ambient_chem)

  return micro


def _micro_step(micro, state, info, opts):
  '''Microphysics step for lagrangian scheme'''
  libopts = lgrngn.opts_t()
  libopts.cond = True
  libopts.coal = False
  libopts.adve = False
  libopts.sedi = False
  libopts.ice_nucl = opts["ice_nucl"]
  libopts.depo = opts["depo"]

  # chemical options
  if micro.opts_init.chem_switch:
    # chem processes: dissolving, dissociation, reactions
    libopts.chem_dsl = opts["chem_dsl"]
    libopts.chem_dsc = opts["chem_dsc"]
    libopts.chem_rct = opts["chem_rct"]

  # get trace gases
  ambient_chem = {}
  if micro.opts_init.chem_switch:
    ambient_chem = dict((v, state[k]) for k,v in _Chem_g_id.items())

  # call libcloudphxx microphysics
  micro.step_sync(libopts, state["th_d"], state["r_v"], state["rhod"], ambient_chem=ambient_chem)
  micro.step_async(libopts)

  # update state after microphysics (needed for below update for chemistry)
  _stats(state, info)

  # update in state for aqueous chem (TODO do we still want to have aq chem in state?)
  if micro.opts_init.chem_switch:
    micro.diag_all() # selecting all particles
    for id_str, id_int in _Chem_g_id.items():
      # save changes due to chemistry
      micro.diag_chem(id_int)
      state[id_str.replace('_g', '_a')] = np.frombuffer(micro.outbuf())[0]
  if micro.opts_init.ice_switch:
    micro.diag_ice()
    micro.diag_ice_mix_ratio()
    state["ice_mix_ratio"] = np.frombuffer(micro.outbuf())[0]


  micro.diag_rw_ge_rc()
  mom0 = micro.diag_wet_mom(0)
  mom0 = np.frombuffer(micro.outbuf())[0]
  state["act_m0"] = mom0

  micro.diag_rw_ge_rc()
  mom1 = micro.diag_wet_mom(1)
  mom1 = np.frombuffer(micro.outbuf())[0]
  state["act_m1"] = mom1

  micro.diag_rw_ge_rc()
  mom2 = micro.diag_wet_mom(2)
  mom2 = np.frombuffer(micro.outbuf())[0]
  state["act_m2"] = mom2

  micro.diag_rw_ge_rc()
  mom3 = micro.diag_wet_mom(3)
  mom3 = np.frombuffer(micro.outbuf())[0]
  state["act_m3"] = mom3

  micro.diag_rw_ge_rc()
  mom4 = micro.diag_wet_mom(4)
  mom4 = np.frombuffer(micro.outbuf())[0]
  state["act_m4"] = mom4

  micro.diag_all()
  micro.diag_sd_conc()
  state["sd_conc"] = np.frombuffer(micro.outbuf())[0]