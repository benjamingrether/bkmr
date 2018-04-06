# -*- coding: utf-8 -*-
"""
Created on Mon Mar 26 12:41:33 2018

@author: bmg18
"""

### Check what is going on with ztest

import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from scipy.special import gamma

def rdelta_comp_update(r, delta, lambda_, y, X, beta, sigsq_eps, Vcomps, Z, ztest, data_comps, control_params, rprop_gen2, rprop_logdens1, rprior_logdens, rprior_logdens2, rprop_logdens2, rprop_gen1):
    r_params = control_params['r_params']
    a_p0 = control_params['a_p0']
    b_p0 = control_params['b_p0']
    delta_star = delta
    r_star = r
    
    if all(delta==0):
        move_type = 1
        move_prob = 1
    else:
        move_type = np.random.choice([1,2])
        move_prob = 0.5
    
    if (move_type == 1):
        comp = np.random.choice(ztest)
        r_params = set_r_params(r_prior = control_params['r_prior'], comp = comp, r_params = r_params)
        
        delta_star[comp] = 1.0 - delta[comp]
        if all(delta_star == 0):
            move_prob_star = 1.0
        else:
            move_prob_star = 0.5
        if (delta_star[comp] == 0):
            r_star[comp] = 0.0
        else:
            r_star[comp] = rprop_gen1(r_params = r_params)
        
        if delta[comp]==1:
            deltacomp_f=-1
            r_sel = r[comp]
        else:
            deltacomp_f=1
            r_sel=r_star[comp]
        
        diffpriors = np.log(gamma(np.sum(delta_star) + a_p0[0])) + np.log(gamma(len(ztest) - np.sum(delta_star) + b_p0[0])) - np.log(gamma(np.sum(delta) + a_p0[0])) - np.log(gamma(len(ztest) - np.sum(delta) + b_p0[0])) + deltacomp_f * rprior_logdens(x = r_sel, r_params=r_params)
        
        negdifflogproposal = (-1.0)*np.log(move_prob_star)+np.log(move_prob)*deltacomp_f*rprop_logdens1(x=r_sel, r_params = r_params)
        
    elif (move_type == 2):
        comp = np.random.choice([i for i, x in enumerate(delta_star) if x == 1])
        r_params = set_r_params(r_prior = control_params['r_prior'], comp = comp, r_params = r_params)
        
        r_star[comp] = rprop_gen2(current = r[comp], r_params = r_params)
        
        diffpriors = rprior_logdens(x=r_star[comp], r_params=r_params) - rprior_logdens(x=r[comp], r_params=r_params)
        
        negdifflogproposal = (-1.0)*rprop_logdens2(prop = r_star[comp], current = r[comp], r_params = r_params) + rprop_logdens2(prop = r[comp], current = r_star[comp], r_params = r_params)
    
    lambda_star = lambda_
    
    # M-H step
    return MHstep(r=r, lambda_=lambda_, lambda_star=lambda_star, r_star=r_star, delta=delta, delta_star=delta_star, y=y, X=X, Z=Z, beta=beta, sigsq_eps=sigsq_eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, Vcomps=Vcomps, move_type=move_type, data_comps=data_comps)






