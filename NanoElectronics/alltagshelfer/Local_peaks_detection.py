#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 28 16:33:30 2020

@author: filipk
"""

import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
from scipy.signal import savgol_filter
from scipy.signal import find_peaks

from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{url}\usepackage{sansmath}\sansmath\usepackage[russian]{babel}') 


df=pd.read_csv('/home/filipk/Desktop/Time_evolution_cor_final/Co/1/Len96at[84-179].csv')
df=pd.read_csv('/home/filipk/Desktop/Time_evolution_cor_final/Co/1/Len278at[259-536].csv')
df=pd.read_csv('/home/filipk/Desktop/Time_evolution_cor_final/Co/0/Len29at[17-45].csv')
df=pd.read_csv('/home/filipk/Desktop/Time_evolution_cor_final/Co/1/Len96at[84-179].csv')
dist=list(df['distance'])
coupling=list(df['coupling'])

#plt.scatter(dist,coupling)
dist=dist[:-2]
coupling=coupling[:-2]

"This part interpolates the data so it is equaly distanced for smoothing, and then uses savitzky-golay filtering"
xnew = np.linspace(dist[0], dist[-1], num=((len(coupling)//10)+1)*10, endpoint=True)
f=interp1d(dist, coupling, kind='linear')
yhat = savgol_filter(f(xnew), 15, 3)
f=interp1d(xnew, yhat, kind='linear')
yhat = savgol_filter(f(xnew), 15, 3)
plt.scatter(dist, coupling,s=2)
plt.plot(xnew,yhat,color='orange')
peaks, _ = find_peaks(yhat, height=0,distance=6)
plt.plot(xnew[peaks], yhat[peaks], "x",color='red')


"This part plots vertical lines from aprox bottom of the plot to the detected maxima"
bottom, top = plt.ylim()
plt.vlines(xnew[peaks],ymin=bottom,ymax=yhat[peaks],linestyle='dashed')
#plt.arrow(xnew[peaks[1]],bottom,xnew[peaks[1]]+1,bottom)

"The next part marks the arrows <-> and puts them from coordinate xy to coordinate xytext, "
"then it plots the distance in the middle of arrow with center aligment"
for i in range(0,len(peaks)-1):
    plt.annotate(s='', xy=(xnew[peaks[i]],bottom), xytext=(xnew[peaks[i+1]],bottom), arrowprops=dict(arrowstyle='<->'))
    plt.text((xnew[peaks[i+1]]+xnew[peaks[i]])/2,bottom*1.1,s="{0:.2f}".format(xnew[peaks[i+1]]-xnew[peaks[i]]),horizontalalignment='center')

"This part puts the xlabel and ylabel in latex format"
plt.xlabel(r'Distance[\AA]',size=14)
plt.ylabel(r'$\mathsf{\Gamma}\, {[\mathsf{eV}]}$',size=14)
plt.savefig('/home/filipk/Desktop/Evol_max.png',dpi=400)