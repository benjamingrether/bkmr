# -*- coding: utf-8 -*-
"""
Created on Wed Apr 04 17:32:11 2018

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

# Run setting
kmbayes(y=pd.DataFrame(np.random.random(10)),Z=pd.DataFrame(np.random.random((10,3))),
        X=pd.DataFrame(np.random.random((10,2))),id=[1,1,2,3,4,1,1,3,3,4],varsel=True)

# R
# kmbayes(y=rnorm(10),Z=matrix(rnorm(30),ncol=3),X=matrix(rnorm(20),ncol=2),
# id=c(1,1,2,3,4,1,1,3,3,4),iter=1000,varsel=FALSE)