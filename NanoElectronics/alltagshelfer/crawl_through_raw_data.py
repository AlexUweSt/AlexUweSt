# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 14:20:52 2020

@author: strobe46
"""

import numpy as np
import os
import matplotlib.pyplot as plt
#%%    
directory = 'C:/Users/Strobel/Nextcloud/Evaluation_data/' + \
    'IETS_measurements_revised/2020_measurements/S3_2020/C60/'+ \
        'E11/13_11/60/'

#%
#for fi in os.listdir(mainpath):
a = []
for file in os.listdir(directory):
    a.append(file)
b = [k for k in a if "IETS" in k] # or "upsweep"or "downsweep"
#print(b)
j = [os.path.join(directory,x) for x in b]
#%
NDRs = []
Lias = []
Voltages = []
OCs = []
for idx, file in enumerate(j):
    if os.path.getsize(file) < 1.2E3:
        pass
    else:
        print(idx)
        MEAS = np.genfromtxt(file,skip_header=8, delimiter = '\t',\
                          usecols = (0,1,2,3,4,5,6,7,8), dtype=np.str)
        MEAS = np.char.replace(MEAS, ',', '.')
        MEAS = np.char.replace(MEAS, '\'', '')
        MEAS = np.char.replace(MEAS, 'b', '').astype(np.float64)
        meas_t = np.transpose(MEAS)
        LIA1 = meas_t[6]
        Header=np.genfromtxt(file,dtype='str',max_rows=5,\
                     usecols = (0,1,2,3,4,5,6,7,8))
        Header = np.char.replace(Header, ',', '.')
        if min(LIA1) < -1E-2:
            OCnum = int(Header[3][1])
            Bias = meas_t[8]
            Lias.append(LIA1)
            Voltages.append(Bias)
            NDRs.append((file[-3:],OCnum, min(LIA1), file))
#%%

number = int(NDRs[0][0])

