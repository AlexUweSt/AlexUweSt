#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 10 11:33:38 2021

@author: Strobe46
"""


import matplotlib.pyplot as plt
import numpy as np


def ad_sub_scatter(V,y,fig,ax,plotpos,x_label,y_label):
    """ads a subplot to figure 'fif' at postition 'plotpos'
    plotpos is given as triple 'eg. 111' or [111] with the axis name ax
    """
    ax = fig.add_subplot(plotpos)
    c = ['k','b','m','c']
    #s = [0.8,1,1,1]
    for i in range(len(V)):
        ax.scatter(V[i],y[i],color=c[i], s=5)
    ax.set_xlabel(x_label, fontsize = 16)
    ax.set_ylabel(y_label, fontsize = 16)
    ax.tick_params(axis='both' , labelsize = 12.0)
    return ax

def fig_scatter(x,y,xlab, ylab):
    f_sc = plt.figure(1,figsize=[6.4, 4.8])
    ad_sub_scatter(x,y,f_sc,'ax1',111,xlab, ylab)
    f_sc.tight_layout()
    return f_sc


#%%

leng = [2.46, 4.49, 2.45, 2.70, 4.73]
wid = [273, 280, 110, 59, 76]
pepo = [60, 116, 207, 639, 625]
T = [7,9,6,4, 4]
device= [2,3,5,10, 11]


f = fig_scatter(leng, pepo, 'length [\u00b5m]', 'Peakposition [mV]')


