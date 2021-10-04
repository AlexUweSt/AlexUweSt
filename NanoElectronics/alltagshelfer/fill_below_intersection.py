#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  8 11:31:50 2021

@author: filipk
"""







#%%

import matplotlib.pyplot as plt
import numpy as np

def boltzman(x, xmid, tau):
    """
    evaluate the boltzman function with midpoint xmid and time constant tau
    over x
    """
    return 1. / (1. + np.exp(-(x-xmid)/tau))

x = np.arange(-6, 6, .01)
S = boltzman(x, 0, 1)
Z = 1-boltzman(x, 0.5, 1)
plt.plot(x, S, x, Z, color='red', lw=2)
plt.show()

#%%

import matplotlib.pyplot as plt
import numpy as np

def boltzman(x, xmid, tau):
    """
    evaluate the boltzman function with midpoint xmid and time constant tau
    over x
    """
    return 1. / (1. + np.exp(-(x-xmid)/tau))

def fill_below_intersection(x, S, Z):
    """
    fill the region below the intersection of S and Z
    source: https://scipy-cookbook.readthedocs.io/items/Matplotlib_SigmoidalFunctions.html
    """
    #find the intersection point
    ind = int(np.nonzero(np.absolute(S-Z)==min(np.absolute(S-Z)))[0])
    # compute a new curve which we will fill below
    Y = np.zeros(S.shape, dtype=np.float_)
    Y[:ind] = S[:ind]  # Y is S up to the intersection
    Y[ind:] = Z[ind:]  # and Z beyond it
    plt.fill(x, Y, facecolor='blue', alpha=0.5)

x = np.arange(-6, 6, .01)
S = boltzman(x, 0, 1)
Z = 1-boltzman(x, 0.5, 1)
f = plt.figure()
ax = f.add_subplot(111)
ax.plot(x, S, x, Z, color='red', lw=2)
fill_below_intersection(x, S, Z)
plt.show()
