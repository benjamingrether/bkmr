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


def set_r_MH_functions(r_prior):
    if r_prior == 'invunif':
        def rprior_logdens(x, r_params):
            r_a = r_params['r_a']
            r_b = r_params['r_b']
            if (1/r_b <= x) and (x <= 1/r_a):
                return -2*np.log(x) - np.log(r_b - r_a)
            else:
                return np.log(0)
        def rprop_gen1(r_params):
            r_a = r_params['r_a']
            r_b = r_params['r_b']
            return 1/np.random.uniform(r_a,r_b)      
        def rprop_logdens1(x, r_params):
            r_a = r_params['r_a']
            r_b = r_params['r_b']
            if (1/r_b <= x) and (x <= 1/r_a):
                return -2*np.log(x) - np.log(r_b - r_a)
            else:
                return np.log(0)     
        def rprop_gen2(current, r_params):
            r_a = r_params['r_a']
            r_b = r_params['r_b']
            r_jump = r_params['r_jump2']
            return stats.truncnorm.rvs(size=1, a=1/r_b, b=1/r_a, loc=current, scale=r_jump)
        def rprop_logdens2(prop, current, r_params):
            r_a = r_params['r_a']
            r_b = r_params['r_b']
            r_jump = r_params['r_jump2']
            return np.log(stats.truncnorm.pdf(x=prop, a=1/r_b, b=1/r_a, loc=current, scale=r_jump))
        def rprop_gen(current, r_params):
            r_a = r_params['r_a']
            r_b = r_params['r_b']
            r_jump = r_params['r_jump']
            return stats.truncnorm.rvs(size=1, a=1/r_b, b=1/r_a, loc=current, scale=r_jump)
        def rprop_logdens(prop, current, r_params):
            r_a = r_params['r_a']
            r_b = r_params['r_b']
            r_jump = r_params['r_jump']
            return np.log(stats.truncnorm.pdf(x=prop, a=1/r_b, b=1/r_a, loc=current, scale=r_jump))
    return {'rprior_logdens': rprior_logdens, 'rprop_gen1': rprop_gen1, 
            'rprop_logdens1': rprop_logdens1, 'rprop_gen2': rprop_gen2,
            'rprop_logdens2': rprop_logdens2, 'rprop_gen': rprop_gen,
            'rprop_logdens': rprop_logdens} 



