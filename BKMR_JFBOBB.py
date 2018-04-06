# -*- coding: utf-8 -*-
"""
Created on Thu Mar 22 16:14:17 2018

@author: bmg18
"""


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import pandas as pd
from scipy import stats
import statsmodels.api as sm
import datetime


def sigmoid(x,location,scale):
  return 1 / (1 + math.exp(-(x-location)/scale))


def SimData(n=50, M=4, sigsq_true=0.5, beta_true=2):
    cov = np.array([[ 0.72, 0.65, 0.45, 0.48],
                    [ 0.65, 0.78, 0.48, 0.55],
                    [ 0.45, 0.48, 0.56, 0.43],
                    [ 0.48, 0.55, 0.43, 0.71]])
    
    np.random.seed(111)
    Z = np.random.multivariate_normal(np.zeros(M), cov,size=n)
    z1 = np.linspace(min(Z[:,1]), max(Z[:,1]), 20)
    z2 = np.linspace(min(Z[:,2]), max(Z[:,2]), 20)
    z = np.vstack((z1,z2)).T
    X = 3*np.cos(Z[:,0])+2*np.random.normal(size=n)   
    eps = np.random.normal(scale=sigsq_true,size=n)
    h = np.zeros(n)
    for i in range(n):
        h[i] = 4*sigmoid(0.25*(Z[i,0]+Z[i,1]+0.5*Z[i,0]*Z[i,1]),0,0.3)
    h_grid = np.zeros((20,20))
    for i in range(20):
        for j in range(20):
            h_grid[i,j]=4*sigmoid(0.25*(z[i,0]+z[j,1]+0.5*z[i,0]*z[j,1]),0,0.3)
    
    mu = beta_true * X + h
    y = mu + eps
        
    return [n, M, sigsq_true, beta_true, Z, z, z1, z2, h, h_grid, X, y]



#dat = SimData(n=50)
#z1 = dat[6]
#z2 = dat[7]
#h_grid = dat[9]
#fig = plt.figure()
#ax = plt.axes(projection='3d')
#ax.plot_surface(z1,z2,h_grid)

#Z = dat[4]
#X = dat[10]
#y = dat[11]

# y: a vector of outcome data of length \code{n}.
# Z: an \code{n}-by-\code{M} matrix of predictor variables to be included in the \code{h} function. Each row represents an observation and each column represents an predictor.
# X: an \code{n}-by-\code{K} matrix of covariate data where each row represents an observation and each column represents a covariate. Should not contain an intercept column.
# iter: number of iterations to run the sampler
# id: optional vector (of length \code{n}) of grouping factors for fitting a model with a random intercept. If NULL then no random intercept will be included.
# verbose: TRUE or FALSE: flag indicating whether to print intermediate diagnostic information during the model fitting.
# Znew: optional matrix of new predictor values at which to predict \code{h}, where each row represents a new observation. This will slow down the model fitting, and can be done as a post-processing step using \code{\link{SamplePred}}
# starting_values: list of starting values for each parameter. If not specified default values will be chosen.
# control_params: list of parameters specifying the prior distributions and tuning parameters for the MCMC algorithm. If not specified default values will be chosen.
# varsel: TRUE or FALSE: indicator for whether to conduct variable selection on the Z variables in \code{h}
# groups: optional vector (of length \code{M}) of group indictors for fitting hierarchical variable selection if varsel=TRUE. If varsel=TRUE without group specification, component-wise variable selections will be performed. Has to be pd.Series type # groups=pd.Series([0,0,1,2,1])
# rmethod: for those predictors being forced into the \code{h} function, the method for sampling the \code{r[m]} values. Takes the value of 'varying' to allow separate \code{r[m]} for each predictor; 'equal' to force the same \code{r[m]} for each predictor; or 'fixed' to fix the \code{r[m]} to their starting values
# est.h: TRUE or FALSE: indicator for whether to sample from the posterior distribution of the subject-specific effects h_i within the main sampler. This will slow down the model fitting.
# return: an object of class "bkmrfit", which has the associated methods:
def kmbayes(y, Z, X = None, iter = 11, id = None, verbose = True, Znew = None, starting_values = None, control_params = None, varsel = False, groups = None, rmethod = "equal", est_h = False):
    
    # Argument check: check vector/matrix sizes 
    if groups is not None:
        if not len(groups)==Z.shape[1]:
            print("Error: the length of the groups vector should be equal to the number of columns of Z")
            return None
    
    nsamp = iter
    
    if X is None:
        X = pd.DataFrame(np.zeros((y.shape[0],1)))
        missingX = True
    else:
        missingX = False
    
    # if id is not empty, build random intercept model
    if id is None: # No id, no random intercept
        randint = False
        nlambda = 1
        crossTT = 0
    else: # there is id, include random intercept
        randint = True
        id = pd.DataFrame(id)
        nid = len(np.unique(id))
        nlambda = 2
        if not (nid == np.max(id))[0]:
            print("Error: id should be a vector max(id)=number of unique ids")
            return None
        # matrix that multiplies the random intercept
        TT = pd.DataFrame(np.zeros((len(id),nid)))
        for i in range(nid):
            for j in range(len(id)):
                TT[i][j]=int(id[0][j]==i+1)
        crossTT = np.dot(TT,np.transpose(TT))#check order in function
        del TT, nid

    data_comps = {'randint': randint, 'nlambda': nlambda, 'crossTT': crossTT}
    del randint, nlambda, crossTT

    # create empty matrices to store the posterior draws in
    chain = {'h_hat': pd.DataFrame(np.zeros((nsamp,Z.shape[0]))),
             'beta': pd.DataFrame(np.zeros((nsamp,X.shape[1]))),
             'lambda': pd.DataFrame(np.full((nsamp,data_comps['nlambda']),np.nan)),
             'sigsq_eps': pd.DataFrame(np.full(nsamp,np.nan)),
             'r': pd.DataFrame(np.full((nsamp,Z.shape[1]),np.nan)),
             'acc_r': pd.DataFrame(np.zeros((nsamp,Z.shape[1]),dtype=bool)),
             'acc_lambda': pd.DataFrame(np.zeros((nsamp,data_comps['nlambda']),dtype=bool)),
             'delta': pd.DataFrame(np.full((nsamp,Z.shape[1]),1.0))}
    
    if varsel:
        chain['acc_rdelta'] = pd.DataFrame(np.full(nsamp,0)) 
        chain['move_type'] = pd.DataFrame(np.full(nsamp,0))
    
    # components to predict h(Znew)
    if Znew is not None:
        chain['hnew'] = pd.DataFrame(np.zeros((nsamp,Znew.shape[0])))
    
    # components if model selection is being done
    if varsel:
        ztest = range(Z.shape[1])
        rdelta_update = rdelta_comp_update
    else:
        ztest = None
    
    # control parameters
    control_params_default = {'lambda_jump': pd.DataFrame(np.full(data_comps['nlambda'],10.0)),
                              'mu_lambda': pd.DataFrame(np.full(data_comps['nlambda'],10.0)),
                              'sigma_lambda': pd.DataFrame(np.full(data_comps['nlambda'],10.0)),
                              'a_p0': [1.0],
                              'b_p0': [1.0],
                              'r_prior': 'invunif',
                              'a_sigsq': [1e-3],
                              'b_sigsq': [1e-3],
                              'mu_r': [5.0],
                              'sigma_r': [5.0],
                              'r_muprop': [1.0],
                              'r_jump': [0.1],
                              'r_jump1': [2.0],
                              'r_jump2': [0.1],
                              'r_a': [0.0],
                              'r_b': [100.0]}
    if control_params is None:
        control_params = control_params_default
    control_params['r_params'] = {'mu_r': control_params['mu_r'],
                                  'sigma_r': control_params['sigma_r'],
                                  'r_muprop': control_params['r_muprop'],
                                  'r_jump': control_params['r_jump'],
                                  'r_jump1': control_params['r_jump1'],
                                  'r_jump2': control_params['r_jump2'],
                                  'r_a': control_params['r_a'],
                                  'r_b': control_params['r_b']}
    
    # components if grouped model selection is being done
    if groups is not None:
        if not varsel:
            print("Error: if doing grouped variable selection, must set varsel = True")
            return None
        rdelta_update = rdelta_group_update
        control_params['group_params'] = {'groups': groups,
                                          'sel_groups': map(lambda x: min(np.where(groups==x)[0]),np.unique(groups)),
                                          'neach_group': map(lambda x: np.sum(groups==x)[0],np.unique(groups))}
    
    # specify functions for doing the Metropolis-Hastings steps to update r
    rfn = set_r_MH_functions(r_prior = control_params['r_prior'])
    rprior_logdens = rfn['rprior_logdens']
    rprop_logdens = rfn['rprop_logdens']
    rprop_logdens1 = rfn['rprop_logdens1']
    rprop_logdens2 = rfn['rprop_logdens2']
    rprop_gen = rfn['rprop_gen']
    rprop_gen1 = rfn['rprop_gen1']
    rprop_gen2 = rfn['rprop_gen2']
    del rfn
    
    # initial values
    starting_values0 = {'h_hat': 1, 'beta': None, 'sigsq_eps': None, 'r': 1, 'lambda': 10, 'delta': 1}
    if starting_values is None:
        starting_values = starting_values0
    if starting_values['beta'] is None or starting_values['sigsq_eps'] is None:
        lmfit0 = sm.OLS(endog=y,exog=sm.add_constant(np.concatenate((X,Z),axis=1))).fit()
        if starting_values['beta'] is None:
            coefX = np.array(lmfit0.params[1:X.shape[1]+1])
            starting_values['beta'] = coefX
        if starting_values['sigsq_eps'] is None:
            starting_values['sigsq_eps'] = lmfit0.mse_resid
    
    # initialise chain
    chain['h_hat'].iloc[0] = starting_values['h_hat']
    chain['beta'].iloc[0] = starting_values['beta']
    chain['lambda'].iloc[0] = starting_values['lambda']
    chain['sigsq_eps'].iloc[0] = starting_values['sigsq_eps']
    chain['r'].iloc[0] = starting_values['r']
    if varsel:
        chain['delta'].iloc[0] = starting_values['delta']
    if groups is not None: # set delta equal 1 for first representative of group; for the others 0
        starting_values['delta'] = np.zeros(len(groups))
        starting_values['delta'][map(lambda x: min(np.where(groups==x)[0]),np.unique(groups))] = 1
        if varsel:
            chain['delta'].iloc[0] = starting_values['delta']
            chain['r'].iloc[0] = chain['delta'].iloc[0]*starting_values['r']
    chain['est_h'] = est_h
    
    # components
    Vcomps = makeVcomps(r = chain['r'].iloc[0], lambda_ = chain['lambda'].iloc[0], Z = Z, data_comps = data_comps)
    
    # start sampling
    chain['time1'] = datetime.datetime.now()
    for s in range(1,nsamp):
        
        # generate posterior samples from marginalized distribution P(beta, sigsq.eps, lambda, r | y)
        
        # beta
        if not missingX:
            chain['beta'].iloc[s] = beta_update(X = X, Vinv = Vcomps['Vinv'], y = y, sigsq_eps = chain['sigsq_eps'].iloc[s-1]).T.values
        
        # sigma_eps^2
        chain['sigsq_eps'].iloc[s] = sigsq_eps_update(y = y, X = X, beta = chain['beta'].iloc[s], Vinv = Vcomps['Vinv'], a_eps = control_params['a_sigsq'], b_eps = control_params['b_sigsq'])
        
        # lambda
        lambdaSim = chain['lambda'].iloc[s-1]
        for comp in range(0,data_comps['nlambda']):
            varcomps = lambda_update(r = chain['r'].iloc[s-1], delta = chain['delta'].iloc[s-1], lambda_ = lambdaSim, whichcomp = comp, y = y, X = X, Z = Z, beta = chain['beta'].iloc[s], sigsq_eps = chain['sigsq_eps'].iloc[s], Vcomps = Vcomps, data_comps = data_comps, control_params = control_params)
            lambdaSim = varcomps['lambda']
            if varcomps['acc']:
                Vcomps = varcomps['Vcomps']
                chain['acc_lambda'].iloc[s][comp] = varcomps['acc']
        chain['lambda'].iloc[s] = lambdaSim.T.values
        
        # r
        rSim = chain['r'].iloc[s-1]
        comp = ztest
        if varsel: # check if this condition is consistent with comp <- which(!1:ncol(Z) %in% ztest), if (length(comp) != 0)
            if (rmethod=="equal"): # common r for those variables not being selected
                varcomps = r_update(r = rSim, whichcomp = comp, delta = chain['delta'].iloc[s-1], lambda_ = chain['lambda'].iloc[s], y = y, X = X, beta = chain['beta'].iloc[s], sigsq_eps = chain['sigsq_eps'].iloc[s], Vcomps = Vcomps, Z = Z, data_comps = data_comps, control_params = control_params, rprior_logdens = rprior_logdens, rprop_gen1 = rprop_gen1, rprop_logdens1 = rprop_logdens1, rprop_gen2 = rprop_gen2, rprop_logdens2 = rprop_logdens2, rprop_gen = rprop_gen, rprop_logdens = rprop_logdens)
                rSim = varcomps['r']
                if varcomps['acc']:
                    Vcomps = varcomps['Vcomps']
                    chain['acc_r'].iloc[s][comp] = varcomps['acc']
            elif (rmethod=="varying"): # allow a different r_m ### CHECK IF ACTUALLY WORKS
                for whichr in comp:
                    varcomps = r_update(r = rSim, whichcomp = whichr, delta = chain['delta'].iloc[s-1], lambda_ = chain['lambda'].iloc[s], y = y, X = X, beta = chain['beta'].iloc[s], sigsq_eps = chain['sigsq_eps'].iloc[s], Vcomps = Vcomps, Z = Z, data_comps = data_comps, control_params = control_params, rprior_logdens = rprior_logdens, rprop_gen1 = rprop_gen1, rprop_logdens1 = rprop_logdens1, rprop_gen2 = rprop_gen2, rprop_logdens2 = rprop_logdens2, rprop_gen = rprop_gen, rprop_logdens = rprop_logdens)
                    rSim = varcomps['r']
                    if varcomps['acc']:
                        Vcomps = varcomps['Vcomps']
                        chain['acc_r'].iloc[s][whichr] = varcomps['acc']
        
        # for those variables being selected: joint posterior of (r,delta)
        if varsel:
            varcomps = rdelta_update(r = rSim, delta = chain['delta'].iloc[s-1], lambda_ = chain['lambda'].iloc[s], y = ycont, X = X, beta = chain['beta'].iloc[s], sigsq_eps = chain['sigsq_eps'].iloc[s], Vcomps = Vcomps, Z = Z, ztest = ztest, data_comps = data_comps, control_params = control_params, rprior_logdens = rprior_logdens, rprop_gen1 = rprop_gen1, rprop_logdens1 = rprop_logdens1, rprop_gen2 = rprop_gen2, rprop_logdens2 = rprop_logdens2, rprop_gen = rprop_gen, rprop_logdens = rprop_logdens)
        
        
        
        
        #### Construction side
        if (s == 2):
            print("over")
            return None
    
    
    #print(chain['sigsq_eps'])
    return chain['acc_lambda']#Vcomps, chain['time1']#, chain['delta']


# kmbayes(y=pd.DataFrame(np.random.random(10)),Z=pd.DataFrame(np.random.random((10,3))),X=pd.DataFrame(np.random.random((10,2))),id=[1,1,2,3,4,1,1,3],varsel=False)
