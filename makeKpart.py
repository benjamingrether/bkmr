# -*- coding: utf-8 -*-
"""
Created on Wed Mar 28 10:34:57 2018

@author: bmg18
"""

import numpy as np
import matplotlib.pyplot as plt
import math
import pandas as pd
from scipy import stats
from scipy.spatial.distance import cdist

def makeKpart(r, Z1, Z2 = None):
    Z1r = Z1.apply(lambda row: row*np.sqrt(r.values.flatten()),axis=1)
    if Z2 is None:
        Z2r = Z1r
    else:
        Z2r = Z2.apply(lambda row: row*np.sqrt(r),axis=1)
    Kpart = cdist(Z1r,Z2r)**2
    return Kpart

