# -*- coding: utf-8 -*-
"""
Created on Wed Apr 04 18:07:22 2018

@author: bmg18
"""

import numpy as np
import pandas as pd

def beta_update(X, Vinv, y, sigsq_eps):
    XVinv = np.dot(np.transpose(X),Vinv)#check order in function
    cholXV = np.linalg.cholesky(np.matmul(XVinv,X))
    Vbeta = np.matmul(np.linalg.inv(cholXV).T,np.linalg.inv(cholXV))
    cholVbeta = np.linalg.cholesky(Vbeta)
    betahat = np.matmul(np.matmul(Vbeta, XVinv), y)
    n01 = np.random.normal(size=X.shape[1])
    return betahat + pd.DataFrame(np.dot(np.transpose(np.sqrt(sigsq_eps)[0]*pd.DataFrame(cholVbeta)),n01))


