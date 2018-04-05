# -*- coding: utf-8 -*-
"""
Created on Thu Apr 05 12:58:08 2018

@author: bmg18
"""


import numpy as np
import pandas as pd
from scipy.stats import gamma

def lambda_update(r, delta, lambda_, whichcomp, y, X, Z, beta, sigsq_eps, Vcomps, data_comps, control_params):
    lambda_jump = control_params['lambda_jump'][0][whichcomp]
    mu_lambda = control_params['mu_lambda'][0][whichcomp]
    sigma_lambda = control_params['sigma_lambda'][0][whichcomp]
    lambdacomp = lambda_[whichcomp]
    
    # generate a proposal
    rate = lambdacomp/(lambda_jump**2)
    lambdacomp_star = np.random.gamma(shape = (lambdacomp/lambda_jump)**2, scale = 1/rate, size = 1)
    r_star = r
    delta_star = delta
    move_type = None
    
    # part of M-H ratio that depends on the proposal distribution
    rate1 = lambdacomp/(lambda_jump**2)
    rate2 = lambdacomp_star/(lambda_jump**2)
    negdifflogproposal = (-1)*np.log(gamma.pdf(lambdacomp_star, a = (lambdacomp/lambda_jump)**2, scale = 1/rate1)) + np.log(gamma.pdf(lambdacomp, a=(lambdacomp_star/lambda_jump)**2, scale=1/rate2)) # check if a is shape
    
    # prior distribution
    rate1 = mu_lambda/(sigma_lambda**2)
    diffpriors = np.log(gamma.pdf(lambdacomp_star, a = (mu_lambda/sigma_lambda)**2, scale = 1/rate1)) - np.log(gamma.pdf(lambdacomp, a = (mu_lambda/sigma_lambda)**2, scale = 1/rate1))
    
    lambda_star = lambda_
    lambda_star[whichcomp] = lambdacomp_star
    
    # M-H step
    return MHstep(r=r, lambda_=lambda_, lambda_star=lambda_star, r_star=r_star, delta=delta, delta_star=delta_star, y=y, X=X, Z=Z, beta=beta, sigsq_eps=sigsq_eps, diffpriors=diffpriors, negdifflogproposal=negdifflogproposal, Vcomps=Vcomps, move_type=move_type, data_comps=data_comps)


