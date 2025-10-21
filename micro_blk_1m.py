#!/usr/bin/env python
import numpy as np
from libcloudphxx import blk_1m
from micro_blk_1m_common import _opts_init_blk_1m_common
from parcel_common import _stats
        
def _opts_init_blk_1m(opts):
    """Initialize options for bulk microphysics scheme"""
    opts_init = blk_1m.opts_t()
    opts_init = _opts_init_blk_1m_common(opts)

    # no ice proceesses:
    opts_init.homA1 = False
    opts_init.homA2 = False
    opts_init.hetA = False
    opts_init.hetB = False 
    opts_init.depA = False 
    opts_init.depB = False
    opts_init.rimA = False
    opts_init.rimB = False 
    opts_init.melA = False
    opts_init.melB = False
    
    return opts_init


def _micro_step_blk_1m(micro_opts, state, info, opts):
  """Microphysics step for bulk scheme without ice processes"""
  # get state variables as numpy arrays
  p = np.asarray(state["p"])
  dot_th_d = np.zeros_like(state["th_d"])
  dot_rv = np.zeros_like(state["r_v"])
  dot_rc = np.zeros_like(state["rc"])
  dot_rr = np.zeros_like(state["rr"])

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
    
  blk_1m.rhs_cellwise_revap(
    micro_opts,
    dot_th_d,
    dot_rv,
    dot_rc,
    dot_rr,
    state["rhod"],
    p,
    state["th_d"],
    state["r_v"],
    state["rc"],
    state["rr"],
    opts["dt"]
  )

  state["th_d"] += dot_th_d * opts["dt"]
  state["r_v"] += dot_rv * opts["dt"]
  state["rc"] += dot_rc * opts["dt"]
  state["rr"] += dot_rr * opts["dt"]
  
  # Update thermodynamic state
  _stats(state, info)

