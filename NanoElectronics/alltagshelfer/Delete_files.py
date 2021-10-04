# -*- coding: utf-8 -*-
"""
Created on Wed Feb 19 18:08:13 2020

@author: strobe46
"""
#  Delete certain files in directory tree 
import numpy as np
import pandas as pd
import os
import shutil



directory = '/home/filipk/Documents/Alex/Evaluation_data/IETS_measurements_revised/2019/Oct_Nov_eval_2021'



delal, fina = [], []
tefi = []
#%%
for dirpath, dirnames, filenames in os.walk(os.path.abspath(directory)):
    delal.append([dirpath, dirnames, filenames]) 
for i , f in enumerate(delal):
    for ind , files in enumerate(delal[i][1]):
        tefi.append(files)
        if 'dia' in files:
            os.remove(delal[i][0]+'/'+files)
#%%            
for dirpath, dirnames, filenames in os.walk(os.path.abspath(directory)):
    fina.append(dirnames)
for i , f in enumerate(fina):
    for ind , files in enumerate(fina[i]):
        if 'dia' in files:
            try:
                shutil.rmtree(delal[i][0]+'/'+files)
            except OSError as e:
                print("Error: %s : %s" % (dir_path, e.strerror))
            #os.remove(delal[i][0]+'/'+files)           
            
            
            
            
            