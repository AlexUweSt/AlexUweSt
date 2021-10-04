# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 13:50:06 2020

@author: strobe46
"""

import numpy as np


def bogen(r,s):
    b = 2*r*np.arcsin(s/(2*r))
    return b

def höhe(r,s):
    h = r-0.5*np.sqrt((4*r**2)-s**2)
    return h

def Segfläche(r,s):
    h = höhe(r,s)
    b = bogen(r,s)
    x = r-h
    A1 = (r*b/2)
    A2 = (s*x/2)
    A = A1 - A2
    return A

def kfläche(r):
    A = r**2*np.pi
    return A

def Fläche(r,s1,s2,deltas):
    Ak = kfläche(r)
    Am = Ak - (Segfläche(r,s1)+Segfläche(r,s2))
    der, des1, des2 = deltas
    Ag = kfläche(r+der) - (Segfläche(r+der,s1+des1)+Segfläche(r+der,s2+des2))
    Akl = kfläche(r-der)- (Segfläche(r-der,s1-des1)+Segfläche(r-der,s2-des2))
    de1 = Ag - Am
    de2 = Am - Akl
    return Am, Ak, Akl, Ag, (de1, de2)

def thickn(A,rho,m):
    d = m/(rho*A[0])*1E6
    dmin = m/(rho*A[3])*1E6
    dmax = m/(rho*A[2])*1E6
    dd1 = d- dmin
    dd2 = dmax -d
    if dd1 < dd2:
        dd = dd2
    else:
        dd = dd1
    return d, dmin, dmax, dd


#%% Berechnung
#rho = 0.00731 # g/mm3
rho = 0.00730 # g/mm3

# waver radius = 100 +- 0.5 mm - 2mm clamping
# flat 1 = 32.5 +- 2.5 mm
# flat 2 = 18.00 +-2 mm
A = Fläche(49, 32.5, 18, (0.25, 2.5, 2))

#%% d calculation

# Wafer 3 :
ml, mg = 6.9363, 6.9682  #leergewicht, gesamtgewicht
delm = round(mg-ml, 5)

# Mass difference Wafer 1 = 0.0280 g
# Mass difference Wafer 2 = 0.0657 g
# Mass difference Wafer 3 = 0.0319 g

thickn = thickn(A, rho, delm)
print(thickn)
## results:

# d Wafer 1 = 513 +- 4 nm
# d Wafer 2 = 1202 +- 9 nm
# d Wafer 3 = 584 +- 4.8 nm


# Parameter verdampfer: rho=7,3 ; Z=0.841 ;














