# -*- coding: utf-8 -*-
"""
Created on Fri Apr 06 11:21:20 2018

@author: bmg18
"""

# if varsel False then comp=1:ncolZ, else comp=None
import numpy as np
import pandas as pd
from scipy.stats import gamma

def r_update(r, whichcomp, delta, lambda_, y, X, beta, sigsq_eps, Vcomps, Z, data_comps, control_params, rprior_logdens, rprop_gen1, rprop_logdens1, rprop_gen2, rprop_logdens2, rprop_gen, rprop_logdens):
    r_params = control_params['r_params']
    rcomp = np.unique(r)
    if (len(rcomp) > 1):
        print("Error: rcomp should only be 1-dimensional")
        return None
    # generate a proposal
    rcomp_star = rprop_gen(current = rcomp, r_params = r_params)
    lambda_star = lambda_
    delta_star = delta
    move_type = None
    
    # part of M-H ratio that depends on the proposal distribution
    negdifflogproposal = (-1.0)*rprop_logdens(prop = rcomp_star, current = rcomp, r_params = r_params) + rprop_logdens(prop = rcomp, current = rcomp_star, r_params = r_params)
    
    # prior distribution
    diffpriors = rprior_logdens(x = rcomp_star, r_params = r_params) - rprior_logdens(x = rcomp, r_params = r_params)
    
    r_star = r
    # r_star = rcomp_star ???
    
    # M-H step
    return MHstep(r=r, lambda_=lambda_, lambda_star=lambda_star, r_star=r_star, delta=delta, delta_star=delta_star, y=y, X=X, Z=Z, beta=beta, sigsq_eps=sigsq_eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, Vcomps=Vcomps, move_type=move_type, data_comps=data_comps)



