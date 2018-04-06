# -*- coding: utf-8 -*-
"""
Created on Fri Apr 06 13:51:03 2018

@author: bmg18
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd

def set_r_params(r_prior, comp, r_params):
    if r_prior == 'gamma':
        if (len(r_params['mu_r']) > 1): r_params['mu_r'] = r_params['mu_r'][comp]
        if (len(r_params['sigma_r']) > 1): r_params['sigma_r'] = r_params['sigma_r'][comp]
        if (len(r_params['r_jump1']) > 1): r_params['r_jump1'] = r_params['r_jump1'][comp]
        if (len(r_params['r_jump2']) > 1): r_params['r_jump2'] = r_params['r_jump2'][comp]
    if (r_prior == 'unif') or (r_prior == 'invunif'):
        if (len(r_params['r_a']) > 1): r_params['r_a'] = r_params['r_a'][comp]
        if (len(r_params['r_b']) > 1): r_params['r_b'] = r_params['r_b'][comp]
        if (len(r_params['r_jump2']) > 1): r_params['r_jump2'] = r_params['r_jump2'][comp]
    return r_params 

