# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 11:25:32 2018

@author: bmg18
"""


### Check what is going on with ztest

import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from scipy.special import gamma

def rdelta_group_update(r, delta, lambda_, y, X, beta, sigsq_eps, Vcomps, Z, ztest, data_comps, control_params, rprop_gen1, rprior_logdens, rprop_logdens1, rprop_gen2, rprop_logdens2):
    r_params = control_params['r_params']
    a_p0 = control_params['a_p0']
    b_p0 = control_params['b_p0']
    groups = control_params['group_params']['groups']
    sel_groups = control_params['group_params']['sel_groups']
    neach_group = pd.DataFrame(control_params['group_params']['neach_group'])
    delta_star = delta
    r_star = r
    
    delta_source = pd.DataFrame(np.zeros(len(np.unique(groups)))) 
    for j in np.unique(groups):
        ind = [i for i, x in enumerate(groups) if x == j-1]
        if any(delta[ind]):
            delta_source[0][j-1] = 1
    delta_source_star = delta_source
    
    # randomly select move type
    if all(delta_source == 0.0):
        move_type = 1
        move_prob = 1
    elif (np.sum((pd.DataFrame(neach_group) > 1.0) == (pd.DataFrame(delta_source) == 1.0))==0)[0]:
        move_type = np.random.choice([1,3])
        move_prob = 0.5
    else:
        move_type = np.random.choice([1,2,3])
        move_prob = 1./3.
    

    if (move_type == 1): # randomly select a source and change its state (e.g., from being in the model to not being in the model)
        
        source = np.random.choice(range(len(delta_source)))
        source_comps = [i for i, x in enumerate(groups) if x == source]
        
        delta_source_star[source] = 1.0 - delta_source[source]
        delta_star[source_comps] = np.random.multinomial(n=1, pvals=np.full(len(source_comps),1.0/len(source_comps)).tolist(), size=int(delta_source_star[0][source]))
        if all(delta_source_star == 0.0):
            move_prob_star = 1
        elif (np.sum((pd.DataFrame(neach_group) > 1.0) == (pd.DataFrame(delta_source_star) == 1.0))==0)[0]:
            move_prob_star = 0.5
        else:
            move_prob_star = 1./3.
        
        # which component got switched
        if delta_source[source]==1:
            inds = [i for i, x in enumerate(delta[source_comps]) if x == 1]
            comps = pd.DataFrame(source_comps)[0][inds]
        else:
            inds = [i for i, x in enumerate(delta_star[source_comps]) if x == 1]
            comps = pd.DataFrame(source_comps)[0][inds]
        r_params = set_r_params(r_prior = control_params['r_prior'], comp = comp, r_params = r_params)
        
        if (delta_star[comp] == 0):
            r_star[comp] = 0.0
        else:
            r_star[comp] = rprop_gen1(r_params = r_params)
        
        if delta_source[source]==1:
            deltacomp_f=1
            r_sel = r[comp]
            diffpriors = np.log(len(sel_groups) - np.sum(delta_source) + b_p0[0]) - np.log(np.sum(delta_source_star + a_p0[0]))
        else:
            deltacomp_f=-1
            r_sel = r_star[comp]
            diffpriors = np.log(np.sum(delta_source) + a_p0[0]) - np.log(len(sel_groups) - np.sum(delta_source_star) + b_p0[0]) +deltacomp_f* np.log(len(source_comps)) +deltacomp_f* rprior_logdens(x = r_sel, r_params = r_params)
        
        negdifflogproposal = (-1.0)*np.log(move_prob_star)+np.log(move_prob)-deltacomp_f*np.log(len(source_comps))-rprop_logdens1(x=r_sel, r_params = r_params)
        
        
    elif (move_type == 2): # randomly select a multi-component source that is in the model and change which component is included
        tmp1 = [i for i, x in enumerate(neach_group) if x > 1]
        tmp2 = [i for i, x in enumerate(delta_source) if x == 1]        
        tmp = common_member(tmp1, tmp2)
        if len(tmp)==0:
            print("Check: tmp lenght==0")
            return None
        if len(tmp)==1:
            source=tmp
        else:
            source=np.random.choice(tmp)
        source_comps = [i for i, x in enumerate(groups) if x == source]
        
        ######
        comp = np.random.choice([i for i, x in enumerate(delta_star) if x == 1])
        r_params = set_r_params(r_prior = control_params['r_prior'], comp = comp, r_params = r_params)
        
        r_star[comp] = rprop_gen2(current = r[comp], r_params = r_params)
        
        diffpriors = rprior_logdens(x=r_star[comp], r_params=r_params) - rprior_logdens(x=r[comp], r_params=r_params)
        
        negdifflogproposal = (-1.0)*rprop_logdens2(prop = r_star[comp], current = r[comp], r_params = r_params) + rprop_logdens2(prop = r[comp], current = r_star[comp], r_params = r_params)
        #######
    
    lambda_star = lambda_
    
    # M-H step
    return MHstep(r=r, lambda_=lambda_, lambda_star=lambda_star, r_star=r_star, delta=delta, delta_star=delta_star, y=y, X=X, Z=Z, beta=beta, sigsq_eps=sigsq_eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, Vcomps=Vcomps, move_type=move_type, data_comps=data_comps)


