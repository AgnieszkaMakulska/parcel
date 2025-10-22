#!/usr/bin/env python
def _opts_init_blk_1m_common(opts_init):    
    opts_init.accr = False  # no rain accretion
    opts_init.conv = False  # no autoconversion
    opts_init.sedi = False  # no sedimentation

    opts_init.adj_nwtrph = True
    opts_init.th_dry = True
    opts_init.const_p = False
        
    return opts_init
