# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 15:38:47 2020

@author: strobe46
"""
#from diffpy.srmise.peaks import TerminationRipples

import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as spc
from scipy.special import wofz
import pylab
def G(x, alpha):
    return np.sqrt(np.log(2) / np.pi)/alpha* np.exp(-(x / alpha)**2 * np.log(2))

def V(x, alpha, gamma):
    sigma = alpha / np.sqrt(2 * np.log(2))
    return np.real(wofz((x + 1j*gamma)/sigma/np.sqrt(2)))/sigma/np.sqrt(2*np.pi)

def L(x, gamma):
    return gamma / np.pi / (x**2 + gamma**2)

# =============================================================================
# def L(v,v0):
#     return (1/2np.pi)/
# =============================================================================

def FWHM(T,V,Wintrinsic=10):
    return np.sqrt(((1.7*V)**2)+((5.4*T*spc.k/spc.e)**2)+Wintrinsic**2)#*1000

# =============================================================================
# # SLM Trasmission calc
# gamma = np.arange(0,5,0.05)
# V = np.arange(-1.4,1.4,0.05)
# e0 = 0.8
# Gamma=0.2
# Tp=((4*Gamma*Gamma)/(((V-e0)**2) + ((Gamma + Gamma)**2)))
# Tn=((4*Gamma*Gamma)/(((V+e0)**2) + ((Gamma + Gamma)**2)))
# plt.close('all')
# plt.plot(V,Tp)
# plt.plot(V,Tn)
# 
# T_ga=[]
# plt.close('all')
# for i in range(len(gamma)):
#     
#     T=(4*gamma[i]*gamma[i])/(((V-e0)**2) + ((gamma[i] + gamma[i])**2))
#     T_ga.append(T)
#     plt.plot(V,T)
#     plt.legend('gamma 0 to 5')
# plt.show()
# =============================================================================

# IETS FWHM influecne:

def Temperature(FWHM_,Vm,Wi=10):
    #T = spc.e*(np.sqrt((FHWM**2)-((1.7*Vm)**2)+(Wi**2)))/(5.4*spc.k)  
    e2 = (spc.e)**2
    V = (1.7*Vm)**2
    W = FWHM_**2
    k = (5.4*spc.k)**2
    Wi2 = Wi**2
    T2 = (e2*(V-W+Wi2))/k
    T_ = np.sqrt(abs(T2))
    T = np.round(T_)
    return T

V = [0.004 , 0.006, 0.008, 0.01, 0.015]
T = np.arange(5,35,1)
VI = 0.004
VII=0.008
VIII=0.015
TI = 10 #in Kelvin

Tcal = Temperature(WI,4)


WI=FWHM(5,4)#VI)
WII=FWHM(10,0.008)
WIII=FWHM(10,0.01)
gamma = 0.01+WI
gammaII = 0.01+WII
gammaIII = 0.01+WIII
LI=L(x, gamma)
LII=L(x, gammaII)
LIII=L(x, gammaIII)
x = np.linspace(-0.8,0.8,10000)
pylab.plot(x, LI, ls='--', label='10 Kelvin Amplitude 4mV')
pylab.plot(x, LII, ls=':', label='10 Kelvin Amplitude 8mV')
pylab.plot(x, LIII, label='10 Kelvin Amplitude 15mV')
pylab.xlim(-0.8,0.8)
pylab.legend()
pylab.show()
WdifT=FWHM(20,VI)-FWHM(5,VI)


go 