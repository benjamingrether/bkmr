# -*- coding: utf-8 -*-
"""
Created on Tue Mar 27 11:27:30 2018

@author: bmg18
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from scipy import stats

### Also implement gamma and unif
def set_r_MH_functions(r_prior):
    if r_prior == 'invunif':
        def rprior_logdens(x, r_params):
            r_a = r_params['r_a']
            r_b = r_params['r_b']
            if (r_a[0] == 0):
                if (1.0/r_b[0] <= x) and (x <= float("inf")):
                    return -2.0*np.log(x) - np.log(r_b[0] - r_a[0])
                else:
                    return float("-inf")
            else:
                if (1.0/r_b[0] <= x) and (x <= 1.0/r_a[0]):
                    return -2.0*np.log(x) - np.log(r_b[0] - r_a[0])
                else:
                    return float("-inf")
        def rprop_gen1(r_params):
            r_a = r_params['r_a']
            r_b = r_params['r_b']
            return 1.0/np.random.uniform(r_a[0],r_b[0])      
        def rprop_logdens1(x, r_params):
            r_a = r_params['r_a']
            r_b = r_params['r_b']
            if (r_a[0] == 0):
                if (1.0/r_b[0] <= x) and (x <= float("inf")):
                    return -2.0*np.log(x) - np.log(r_b[0] - r_a[0])
                else:
                    return float("-inf")
            else:
                if (1.0/r_b[0] <= x) and (x <= 1.0/r_a[0]):
                    return -2.0*np.log(x) - np.log(r_b[0] - r_a[0])
                else:
                    return float("-inf")    
        def rprop_gen2(current, r_params):
            r_a = r_params['r_a']
            r_b = r_params['r_b']
            r_jump = r_params['r_jump2']
            if (r_a[0] == 0):
                return stats.truncnorm.rvs(size=1, a=1.0/r_b[0], b=float("inf"), loc=current, scale=r_jump)
            else:
                return stats.truncnorm.rvs(size=1, a=1.0/r_b[0], b=1.0/r_a[0], loc=current, scale=r_jump)
        def rprop_logdens2(prop, current, r_params):
            r_a = r_params['r_a']
            r_b = r_params['r_b']
            r_jump = r_params['r_jump2']
            if (r_a[0] == 0):
                return np.log(stats.truncnorm.pdf(x=prop, a=1.0/r_b[0], b=float("inf"), loc=current, scale=r_jump))
            else:
                return np.log(stats.truncnorm.pdf(x=prop, a=1.0/r_b[0], b=1.0/r_a[0], loc=current, scale=r_jump))
        def rprop_gen(current, r_params):
            r_a = r_params['r_a']
            r_b = r_params['r_b']
            r_jump = r_params['r_jump']
            if (r_a[0] == 0):
                return stats.truncnorm.rvs(size=1, a=1.0/r_b[0], b=float("inf"), loc=current, scale=r_jump)
            else:
                return stats.truncnorm.rvs(size=1, a=1.0/r_b[0], b=1.0/r_a[0], loc=current, scale=r_jump)
        def rprop_logdens(prop, current, r_params):
            r_a = r_params['r_a']
            r_b = r_params['r_b']
            r_jump = r_params['r_jump']
            if (r_a[0] == 0):
                return np.log(stats.truncnorm.pdf(x=prop, a=1.0/r_b[0], b=float("inf"), loc=current, scale=r_jump))
            else:
                return np.log(stats.truncnorm.pdf(x=prop, a=1.0/r_b[0], b=1.0/r_a[0], loc=current, scale=r_jump))
    return {'rprior_logdens': rprior_logdens, 'rprop_gen1': rprop_gen1, 
            'rprop_logdens1': rprop_logdens1, 'rprop_gen2': rprop_gen2,
            'rprop_logdens2': rprop_logdens2, 'rprop_gen': rprop_gen,
            'rprop_logdens': rprop_logdens} 



