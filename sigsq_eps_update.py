# -*- coding: utf-8 -*-
"""
Created on Thu Apr 05 12:10:14 2018

@author: bmg18
"""

import numpy as np
import pandas as pd

def sigsq_eps_update(y, X, beta, Vinv, a_eps, b_eps):
    mu = y - pd.DataFrame(np.matmul(X,beta))
    rate = b_eps + 0.5*np.matmul(np.dot(np.transpose(mu),Vinv),mu)
    prec_y = np.random.gamma(shape = a_eps + X.shape[0]/2, scale = 1/rate[0], size = 1)
    return 1/prec_y


