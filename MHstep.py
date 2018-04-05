# -*- coding: utf-8 -*-
"""
Created on Thu Apr 05 16:19:56 2018

@author: bmg18
"""

import numpy as np
import pandas as pd

def MHstep(r, lambda_, lambda_star, r_star, delta, delta_star, y, X, Z, beta, sigsq_eps, diffpriors, negdifflogproposal, Vcomps, move_type, data_comps):
    # compute log M-H ratio
    Vcomps_star = makeVcomps(r=r_star, lambda_=lambda_star, Z=Z, data_comps=data_comps)
    mu = y - pd.DataFrame(np.matmul(X,beta))
    diffliks = (1.0/2.0) * Vcomps_star['logdetVinv'] - (1.0/2.0) *  Vcomps['logdetVinv'] - (1.0/2.0/sigsq_eps.values) * np.matmul( np.dot(np.transpose(mu), Vcomps_star['Vinv'] - Vcomps['Vinv']), mu)
    logMHratio = diffliks + diffpriors + negdifflogproposal
    logalpha = min(0,logMHratio)
    
    # return value
    acc = False
    if (np.log(np.random.uniform(size=1)) <= logalpha):
        r = r_star
        delta = delta_star
        lambda_ = lambda_star
        Vcomps = Vcomps_star
        acc = True
    
    return {'r': r, 'lambda': lambda_, 'delta': delta, 'acc': acc, 'Vcomps': Vcomps, 'move_type': move_type}

