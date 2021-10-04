# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 09:29:30 2019

@author: strobe46
"""
#import numpy as np

#### for the peakdetection ############
file = open('Peakparameters.txt' , 'w+')
header ='Peakdiff' '\t' 'ACamp' '\t' 'G0' '\t' 'smooth_window' '\t' 'polyorder'  '\t' 'driv_order' '\t' 'delta_smooth' '\t' 'Peak_looka' '\t' 'delta_peak'
file.write('%s\t\r' % header)
file.close()
###################################################################################

import os
# %% rename files

def rename(old, replo, repln):
    for i,j in enumerate(old):
        new = j.replace(replo, repln)
        os.rename(j, new)

#%%
a =[]
root = 'C:/Users/strobe46/ownCloud/Evaluation_data/IETS_measurements_revised/'

adpath = '2020_measurements/S3_2020/C60/E11/25_11/75/'
directory = os.path.join(root, adpath)

for file in os.listdir(directory):
    a.append(file)
b = [k for k in a if "IETS" in k] # or "upsweep"or "downsweep"
print(b)
c = [k for k in a if "upsweep" in k]
print(c)
d = [k for k in a if "downsweep" in k]
print(d)
e = [k for k in a if "G0" in k]
j = [os.path.join(directory,x) for x in b]   # fullsweep Matrix paths
l = [os.path.join(directory,x) for x in c]    # firstsweep ('upsweep' - Matrix) paths
m = [os.path.join(directory,x) for x in d]   #secondsweep ('downsweep' - Matrix) paths
o = [os.path.join(directory,x) for x in e]   # GO paths

#%%

rename(j, '757-6hz', '187Hz_femt')
rename(l, '757-6hz', '187Hz_femt')
rename(m, '757-6hz', '187Hz_femt')
rename(o, '757-6hz', '187Hz_femt')










