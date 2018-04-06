# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 10:32:44 2018

@author: bmg18
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from scipy import stats

def makeVcomps(r, lambda_, Z, data_comps):
    Kpart = makeKpart(r=r, Z1=Z)
    V = pd.DataFrame(np.identity(Z.shape[0])) + lambda_[0]*pd.DataFrame(np.exp((-1)*Kpart))
    if data_comps['nlambda'] == 2:
        V = V + lambda_[1]*data_comps['crossTT']
    cholV = np.linalg.cholesky(V)
    Vinv = np.matmul(np.linalg.inv(cholV).T,np.linalg.inv(cholV))
    logdetVinv = -2.0*sum(np.log(np.diagonal(cholV)))
    Vcomps = {'Vinv': Vinv, 'logdetVinv': logdetVinv}
    return Vcomps

#makeVcomps(r=pd.DataFrame([1,0,1,1]), lambda_=[1,2], Z=pd.DataFrame(np.identity(4)), data_comps=data_comps)
#makeKpart(r=pd.DataFrame([1,0,1]), Z1=pd.DataFrame(np.identity(3)))
