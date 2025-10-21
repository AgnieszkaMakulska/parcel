#!/usr/bin/env python
import numpy as np
from libcloudphxx import blk_1m
from micro_blk_1m_common import _opts_init_blk_1m_common
from parcel_common import _stats
        
def _opts_init_blk_1m_ice(opts):
    """Initialize options for bulk microphysics scheme with ice processes"""
    opts_init = blk_1m.opts_t()
    opts_init = _opts_init_blk_1m_common(opts_init)
    
    return opts_init

def _micro_step_blk_1m_ice(micro_opts, state, info, opts):
  """Microphysics step for bulk scheme with ice processes"""
  # get state variables as numpy arrays
  p = np.asarray(state["p"])
  dot_th_d = np.zeros_like(state["th_d"])
  dot_rv = np.zeros_like(state["r_v"])
  dot_rc = np.zeros_like(state["rc"])
  dot_rr = np.zeros_like(state["rr"])
  dot_ria = np.zeros_like(state["ria"])
  dot_rib = np.zeros_like(state["rib"])

  blk_1m.adj_cellwise(
    micro_opts,
    state["rhod"],
    p,
    state["th_d"],
    state["r_v"],
    state["rc"],
    state["rr"],
    opts["dt"]
    )
  
  blk_1m.rhs_cellwise_ice(
    micro_opts,
    dot_th_d,
    dot_rv,
    dot_rc,
    dot_rr,
    dot_ria,
    dot_rib,
    state["rhod"],
    p,
    state["th_d"],
    state["r_v"],
    state["rc"],
    state["rr"],
    state["ria"],
    state["rib"],
    opts["dt"]
  )

  state["th_d"] += dot_th_d * opts["dt"]
  state["r_v"] += dot_rv * opts["dt"]
  state["rc"] += dot_rc * opts["dt"]
  state["rr"] += dot_rr * opts["dt"]
  state["ria"] += dot_ria * opts["dt"]
  state["rib"] += dot_rib * opts["dt"]

  # Update thermodynamic state
  _stats(state, info)
