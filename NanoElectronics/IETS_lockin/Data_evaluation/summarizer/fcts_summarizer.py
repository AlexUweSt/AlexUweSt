#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  6 13:06:05 2021

@author: Strobel Alexander

the Functions for the summarizer programm

due: 16.04.2021

"""

import os
#import peakdetect as pedet
#import scipy.constants as spc
#import pickle
#from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.collections import LineCollection
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
#import matplotlib.patches as mpatches
import numpy as np
import json
from matplotlib import cm
from scipy import signal
from lmfit.models import LorentzianModel, QuadraticModel, GaussianModel,\
                                LinearModel

cols = ['r','c','g' ,'b','y','lime', 'cadetblue','m','lightgreen', \
      'dimgrey', 'goldenrod', 'slateblue', 'gold', 'deeppink', \
          'steelblue','k', 'orangered', 'magenta', 'maroon', 'saddlebrown',\
              'salmon', 'tan', 'chartreuse', 'violet', 'teal', 'navy']

def key_dct():  # key_dictionary
    """
    dictionary : longname, dict dicider, keys, label
    """
    dct = {'GI' : ('Integrated_conductance', [0], \
                   ['Integrals', 'didv_whole_range_', 1] , \
                     '$G_{\u222b} \ [G/G_{0}]$'),  
           'G0' : ('G0V', [0], ['G_values', 'Gd0'], '$G_{0V} \ [G/G_{0}]$'), \
               
           'GIn' : ('Integ_neg', [0], ['Integrals', 'didv_neg_arm_', 1], 
                    '$G_{neg\u222b} \ [G/G_{0}]$'), \
               
           'GIp' : ('Integ_pos', [0], ['Integrals', 'didv_pos_arm_', 1], 
                    '$G_{pos\u222b} \ [G/G_{0}]$'), \
               
           'GIr' : ('G_ratio_posneg', [0], ['Integrals', 'didv_ratio_pn_'],  
                    'ratio $G_{pos\u222b} \ / \ G_{neg\u222b}$'), 
           
           'G50' :('G50mV',[0],['G_values','Gd50'],'$G_{50 mV} \ [G/G_{0}]$'),\
               
           'G100' : ('G100mV', [0], ['G_values','Gd100'], 
                     '$G_{100 mV} \ [G/G_{0}]$'), \

           'Gm50' :('Gm50mV',[0],['G_values','Gd_neg_50'], \
                    '$G_{-50 mV} \ [G/G_{0}]$'),\
               
           'Gm100' : ('Gm100mV', [0], ['G_values','Gd_neg_100'], 
                     '$G_{-100 mV} \ [G/G_{0}]$'), \
               
           'Gma' : ('Gmax', [0], ['curve_symmetries', 'dIdV', \
                             'dIdV_max_global', 1], '$G_{max} \ [G/G_{0}]$'), \
               
           'Gmi' : ('Gmin', [0], ['curve_symmetries', 'dIdV', \
                              'dIdV_min_global', 1], '$G_{min} \ [G/G_{0}]$'),\
               
           'Gma_x' : ('Bias_Gmax', [0], ['curve_symmetries', 'dIdV', \
                                 'dIdV_max_global', 0], '$U_{max} \ [mV]$'), \
               
           'Gmi_x' : ('Bias_Gmin', [0], ['curve_symmetries', 'dIdV', 
                                 'dIdV_min_global', 0], '$U_{min} \ [mV]$'),
           
           'Ima' : ('Imax', [0], ['curve_symmetries', 'I','Imax_pos', 1],\
                    '$G_{max} \ [G/G_{0}]$'), \
               
           'Imi' : ('Imin', [0], ['curve_symmetries', 'I', 'Imin_neg', 1],\
                    '$G_{min} \ [G/G_{0}]$'),\
               
           'Ima_x' : ('Bias_Imax', [0], ['curve_symmetries', 'I', \
                                 'Imax_pos', 0], '$U_{Bias} \ [mV]$'), \
               
           'Imi_x' : ('Bias_Imin', [0],['curve_symmetries','I','Imin_neg', 0],\
                      '$U_{Bias} \ [mV]$'),   
           'Mpos' : ('Motorposition', [0], ['Header', 'Motor_Pos'],  
                     'Motorposition [a.u.]'), \
               
           'Gfac' : ('G_factor', [0], ['gen_inf', 'G_scale_fac'], \
                     ' G_scale_factor'), \
               
           'Kormi' : ('Kondo_Range_min', [0], ['Kondo','Kondo_window',0], 
                       'Kondo Range min [mV]'), \

           'Korma' : ('Kondo_Range_max', [0], ['Kondo','Kondo_window',1], 
                       'Kondo Range max [mV]'), \
               
           'TK' : ('Kondo_Temp', [0], ['Kondo', 'Kondo_temp_Fano_norm'], 
                       '$T_{K} \ [K]$'),
           
           'AIndr': ('NDR_Integral_area', [0], ['NDR', 'area'], \
                     '$ Area_{NDR} \  [A]$'), \

           'GIndr': ('NDR_Integral', [0], ['NDR', 'area_rang'], \
                     '$G_{\u222b NDR} \  [G/G_{0}]$'), \

           'ndrcrl': ('NDR_zero_cross_left_y', [0], ['NDR','crossings_x_y',\
                             0, 0], '$U_{B} (I_{0}) \ [mV]$'), \

            'ndrcrr': ('NDR_zero_cross_right_x', [0], ['NDR','crossings_x_y', \
                               0, 1],'$U_{B} (I_{0}) \ [mV]$'), \
                
            'Gndr_r': ('GIntegral_ratio_NDR', [0], ['NDR', 'GI_rat_whole'], \
                       '$Ratio \ G_{\u222b} /G_{NDR}$'),
 
            'Gndr_pr': ('G_posInt_ratio_NDR', [0], ['NDR', 'GI_rat_parm'], \
                       '$Ratio \ G_{pos\u222b}/G{NDR}$'),
               
               ## peak dictionary

           'Gpe_x' : ('G_peaks_x', [2], ['didv_extrema', 0, 0], 
                      '$U_{Bias} \ [mV]$'), 
           
           'Gpe_y' : ('G_peaks_height', [2], ['didv_extrema', 0, 1], 
                      '$ height \ [G/G_{0}]$'), 
           
           'Gdi_x' : ('G_dips_x', [2], ['didv_extrema', 1, 0], 
                      '$U_{Bias} \ [mV]$'), 
           
           'Gdi_y' : ('G_dip_height', [2], ['didv_extrema', 1, 1], 
                      '$ height \ [G/G_{0}]$'),  
           
           'dGpa_x' : ('dG_peaks_antisy_x', [2],['Antisymm_Extrema_dgdv',0,0],\
                      '$U_{Bias} \ [mV]$'),  
               
           'dGpa_y' : ('dG_peaks_antisy_y', [2],['Antisymm_Extrema_dgdv',0,1],\
                      '$dG/dV \ [(G/G_{0})*V^{-1}]$'),     
               
           'dGda_x' : ('dG_dips_antisy_x', [2],['Antisymm_Extrema_dgdv',1,0],\
                      '$U_{Bias} \ [mV]$'),  
               
           'dGda_y' : ('dG_dips_antisy_y', [2],['Antisymm_Extrema_dgdv',1,1],\
                      '$dG/dV \ [(G/G_{0})*V^{-1}]$'),  
               
           'Ipa_x' : ('Int_peaks_antisy_x', [2],['Antisymm_Extrema_iets',0,0],\
                      '$U_{Bias} \ [mV]$'),  
               
           'Ipa_y' : ('Int_peaks_antisy_y', [2],['Antisymm_Extrema_iets',0,1],\
                      '$Intensity \ [V^{-1}]$'),     
               
           'Ida_x' : ('Int_dips_antisy_x', [2],['Antisymm_Extrema_iets',1,0],\
                      '$U_{Bias} \ [mV]$'),  
               
           'Ida_y' : ('Int_dips_antisy_y', [2],['Antisymm_Extrema_iets',1,1],\
                      '$Intensity \ [V^{-1}]$'),                 

           'dGpam_x' : ('dG_peaks_asmeas_x', [2],['as_measured_Extrema_dgdv',\
                      0,0],'$U_{Bias} \ [mV]$'),  

           'dGpam_y' : ('dG_peaks_asmeas_y', [2],['as_measured_Extrema_dgdv',\
                      0,1], '$dG/dV \ [(G/G_{0})*V^{-1}]$'),                 

           'dGdam_x' : ('dG_dips_asmeas_x', [2],['as_measured_Extrema_dgdv',\
                      1,0],'$U_{Bias} \ [mV]$'),  

           'dGdam_y' : ('dG_dips_asmeas_y', [2],['as_measured_Extrema_dgdv',\
                      1,1], '$dG/dV \ [(G/G_{0})*V^{-1}]$'),   
               
           'Ipam_x' : ('Int_peaks_asmeas_x', [2],['as_measured_Extrema_iets', \
                              0,0], '$U_{Bias} \ [mV]$'),  
               
           'Ipam_y' : ('Int_peaks_asmeas_y', [2],['as_measured_Extrema_iets', \
                              0,1], '$Intensity \ [V^{-1}]$'),     
               
           'Idam_x' : ('Int_dips_asmeas_x', [2],['as_measured_Extrema_iets', \
                             1,0], '$U_{Bias} \ [mV]$'),  
               
           'Idam_y' : ('Int_dips_asmeas_y', [2],['as_measured_Extrema_iets', \
                             1,1], '$Intensity \ [V^{-1}]$'),     
               
           'dGpbg_x' : ('dG_peaks_Backg_x', [2],['bkg_subs_Extrema_dgdv',\
                      0,0],'$U_{Bias} \ [mV]$'),  

           'dGpbg_y' : ('dG_peaks_Backg_y', [2],['bkg_subs_Extrema_dgdv',\
                      0,1], '$dG/dV \ [(G/G_{0})*V^{-1}]$'),                 

           'dGdbg_x' : ('dG_dips_Backg_x', [2],['bkg_subs_Extrema_dgdv',\
                      1,0],'$U_{Bias} \ [mV]$'),  

           'dGdbg_y' : ('dG_dips_Backg_y', [2],['bkg_subs_Extrema_dgdv',\
                      1,1], '$dG/dV \ [(G/G_{0})*V^{-1}]$'),   
               
           'Ipbg_x' : ('Int_peaks_Backg_x', [2],['bkg_subs_Extrema_iets', \
                              0,0], '$U_{Bias} \ [mV]$'),  
               
           'Ipbg_y' : ('Int_peaks_Backg_y', [2],['bkg_subs_Extrema_iets', \
                              0,1], '$Intensity \ [V^{-1}]$'),     
               
           'Idbg_x' : ('Int_dips_Backg_x', [2],['bkg_subs_Extrema_iets', \
                             1,0], '$U_{Bias} \ [mV]$'),  
               
           'Idbg_y' : ('Int_dips_Backg_y', [2],['bkg_subs_Extrema_iets', \
                             1,1], '$Intensity \ [V^{-1}]$'), 
               
           'dGpdif_x' : ('dG_peaks_peakdif_x', [2],['peakdifff_Extrema_dgdv',\
                      0,0],'$U_{Bias} \ [mV]$'),  

           'dGpdif_y' : ('dG_peaks_peakdif_y', [2],['peakdifff_Extrema_dgdv',\
                      0,1], '$dG/dV \ [(G/G_{0})*V^{-1}]$'),                 

           'dGddif_x' : ('dG_dips_peakdif_x', [2],['peakdifff_Extrema_dgdv',\
                      1,0],'$U_{Bias} \ [mV]$'),  

           'dGddif_y' : ('dG_dips_peakdif_y', [2],['peakdifff_Extrema_dgdv',\
                      1,1], '$dG/dV \ [(G/G_{0})*V^{-1}]$'),   
               
           'Ipdif_x' : ('Int_peaks_peakdif_x', [2],\
            ['peakdiffference_Extrema_iets',0,0], '$U_{Bias} \ [mV]$'),  
               
           'Ipdif_y' : ('Int_peaks_peakdif_y', [2], \
            ['peakdiffference_Extrema_iets', 0,1], '$Intensity \ [V^{-1}]$'),     
               
           'Iddif_x' : ('Int_dips_peakdif_x', [2], \
            ['peakdiffference_Extrema_iets', 1,0], '$U_{Bias} \ [mV]$'),  
               
           'Iddif_y' : ('Int_dips_peakdif_y', [2], \
            ['peakdiffference_Extrema_iets', 1,1], '$Intensity \ [V^{-1}]$'), 
               
           
           # curve dictionary
           
           'Ub' : ('Bias', [1], ['shifted_spectra', 0], '$U_{Bias} \ [mV]$'), \
           'I' : ('Current', [1], ['shifted_spectra', 1], \
                          '$Current \  [\u00b5A]$'), \
           'G' : ('Differential Conductance', [1], ['shifted_spectra', 2], \
                           '$dI/dV  \ [G/G_{0}]$'),
               
           'dG' : ('sec_derivative', [1], ['shifted_spectra', 3], \
                           '$dG/dV \ [(G/G_{0})*V^{-1}]$'),
               
           'iets' : ('iets_spectrum', [1], ['shifted_spectra', 4], \
                           '$Intensity \ [V^{-1}]$')
               
           } 
    return dct

#===================================

def fp_exists(directory):
    if os.path.exists(directory) == True:
        print('\t valid path')
    else:
        print('\n \t The path does not exist check scorrect spelling ')
        

def load_eval(fp):
    """
    fp: file path string

    """
    dir_ar = fp + '/curve_attributes'
    dir_ex = fp + '/extrema_analysed'
    dir_cu = fp + '/curves_analysed'
    a, b, c = [], [], []
    for file in os.listdir(dir_ar):
        a.append(file) # curve attributes
    for file in os.listdir(dir_cu):
        b.append(file) # curves
    for file in os.listdir(dir_ex):
        c.append(file) # Extrema
    j = [os.path.join(dir_ar,x) for x in a] # curve attributes
    jj = [os.path.join(dir_cu,x) for x in b] # curves
    jjj = [os.path.join(dir_ex,x) for x in c] # extrema
    return j, jj ,jjj


def filt_curv_eval(fp, curves):
    attr, curv, extr  = [], [], []
    ja, jc, je = load_eval(fp)
    if curves == None:
        curves = np.arrange(0, 1E4,1).tolist()
    for i,n in enumerate(ja):
        fnum = int(n[-11:-8])
        #num.append(fnum)
        if fnum in curves:
            with open(n, 'r') as file:
                at=file.read()
            att = json.loads(at)#,allow_pickle=True)
            attr.append(att)
    for i,n in enumerate(jc):
        fnum = int(n[-11:-8])
        if fnum in curves:
            with open(n, 'r') as file:
                cu=file.read()
            cur = json.loads(cu)#,allow_pickle=True)
            curv.append(cur)
    for i,n in enumerate(je):
        fnum = int(n[-11:-8])
        if fnum in curves:
            with open(n, 'r') as file:
                ex=file.read()
            ext = json.loads(ex)#,allow_pickle=True)
            extr.append(ext)
    return attr, curv, extr 

def fbk(dlst, key, fdi, inex): # filter by key
    """
    fbk: filter by key. Filters out the measurement exhibiting a certain 
    feature in a dictionary (key in dictionary ) this feature  is given by the 
    'key' argument. That can be 'Coulomb', 'NDr' or something else we are
    looking for. If we need all curve choose a 'key' argument that is there 
    in every evaluation like 'Integrals', 'Header', etc.. 
    
    Parameters
    ----------
    dlst: list with dictionaries
    key : string
        key to filter dictionaries 2 key if double is 'yes'
    fdi : int
        the dictionary to look for filter pararamter
    inex : string
        include ('in', 'or'), exclude ('exc') second key, only if len(key)==
        include: 1.  key1 and key2 must be in dictionarie
                 2. key1 or key2 must be in dict
        exclude: key1 andnot key2 in dictionary

    Returns
    -------
    list
        list of filtered dictionaries
    nams : TYPE
        DESCRIPTION.

    """
    if fdi == 0:
        di = dlst[0]
    elif fdi == 1:
        di = dlst[1]
    elif fdi == 2:
        di = dlst[2]
    at, cu, ex = dlst
    atf, nam, cuf, exf = [], [], [], [] #atributes filtered, names, curves, filtered, extr filt
    #num = []
    if len(key) == 1:
        for k, l in enumerate(di): 
            if key[0] in l.keys():
                #atf.append(l)
                nam.append(l['name'])
    elif len(key) == 2:
        if inex == 'in':
            for k, l in enumerate(di): 
                if key[0] in l.keys() and key[1] in l.keys():
                    #atf.append(l)
                    nam.append(l['name'])            
        if inex == 'exc':
            for k, l in enumerate(di): 
                if key[0] in l.keys() and not key[1] in l.keys():
                    #atf.append(l)
                    nam.append(l['name'])  
        if inex == 'or':
            for k, l in enumerate(di): 
                if key[0] in l.keys() or key[1] in l.keys():
                    #atf.append(l)
                    nam.append(l['name'])  
    for k, l in enumerate(at):
        if l['name'] in nam:
            atf.append(l)           
    for k, l in enumerate(ex):
        if l['name'] in nam:
            exf.append(l)
    for k, l in enumerate(cu):
        if l['name'] in nam:
            cuf.append(l)
    nams = []#[int(i[-6:-3]) for i in nam]
    for c,k in enumerate(nam):
        indf = k.find('IETS0')  # fullsweepname
        nams.append(int(k[indf+4:indf+8]))
    if len(nam) == 0:
        print('check key string /n arrays empty')
    return [atf, cuf, exf], nams

def fbv(dlst, key, comp, thresh):
    di = dlst[0]
    at, cu, ex = dlst
    atf, nam, cuf, exf = [], [], [], [] #atributes filtered, names, curves, filtered, extr filt
    vals = propflt(dlst, key)
    for i,j in enumerate(vals):
        if comp == '>':
            if j > thresh:
                nam.append(di[i]['name'])
        if comp == '<':
            if j < thresh:
                nam.append(di[i]['name'])
    for k, l in enumerate(at):
        if l['name'] in nam:
            atf.append(l)           
    for k, l in enumerate(ex):
        if l['name'] in nam:
            exf.append(l)
    for k, l in enumerate(cu):
        if l['name'] in nam:
            cuf.append(l)
    nams = []#[int(i[-6:-3]) for i in nam]
    for c,k in enumerate(nam):
        indf = k.find('IETS0')  # fullsweepname
        nams.append(int(k[indf+4:indf+8]))
    if len(nam) == 0:
        print('check key string /n arrays empty')
    return [atf, cuf, exf], nams
                
def cu_extract(cuf, key):
    cuar, cuarr = [], []
    for n, k in enumerate(cuf):
        if len(key) == 1:
            cuar.append(k[key[0]])
        elif len(key) == 2:
            cuar.append(k[key[0]][key[1]])
        elif len(key) == 3:
            cuar.append(k[key[0]][key[1]][key[2]])
    return cuar

def fbysn(dcts, swnam): #filter by sweep names
    at, cu, ex = [], [], []
    for i,j in enumerate(dcts[0]):
        na = int(j['name'][-6:-3])
        if na in swnam:
            at.append(j)
            cu.append(dcts[1][i])
            ex.append(dcts[2][i])
    return [at,cu,ex]

def curname(curs):
    if len(curs) < 400:
        custring = str(min(curs)) + '-' + str(max(curs))
    else:
        custring = 'all'
    return custring

def curve_arrays(d, da, spec):
    """
    specify the curves to be taken from the dict d.
    
    d: the dictionary exhibiting the curves
    key: the dictionary key of the curve to be extracted
    da: the dictionary exhibiting the G scale factor (dict additional)
    spec: string for if statement for the exact the spectrum to extract 
    (key to further specify which curves)
    """
    #arr = cu_extract(d, key)
    x, y = [], []

    if spec == 'didv' or spec == 'iv' or spec == 'dgdv' or spec == 'iets' :
        arr = cu_extract(d, ['shifted_spectra'])
        for i, n in enumerate(arr):
            x.append(np.array(n[0]))       
        if spec == 'didv':
            for i, k in enumerate(arr):
                fac = da[i]['gen_inf']['G_scale_fac']
                y.append(np.array(k[2])/fac)
        elif spec == 'dgdv':
            for i, k in enumerate(arr):
                y.append(np.array(k[3]))
        elif spec == 'iets':
            for i, k in enumerate(arr):
                y.append(np.array(k[4]))
        elif spec == 'iv':
            for i, k in enumerate(arr):
                y.append(np.array(k[1]))
                
    elif spec == 'Kondo_fitf' or spec == 'Kondo_Fitb':
        arr = cu_extract(d, ['Kondo', 'Fano'])
        for i, n in enumerate(arr):
            x.append(np.array(n['V-ax']))
        if spec == 'Kondo_fitf':
            for i, k in enumerate(arr):
                y.append([np.array(k['fit_1'])])
        elif spec == 'Kondo_fitb':
            for i, k in enumerate(arr):
                y.append([np.array(k['best_fit'])])
                
    elif spec == 'Kondo_fitfn' or spec == 'Kondo_Fitbn':
        arr = cu_extract(d, ['Kondo','Fano_normalized'])
        for i, n in enumerate(arr):
            x.append(np.array(n['V-ax']))
        if spec == 'Kondo_fitfn':
            for i, k in enumerate(arr):
                y.append(np.array(k['fit_1']))
        elif spec == 'Kondo_fitbn':
            for i, k in enumerate(arr):
                y.append(np.array(k['best_fit']))
    elif spec == 'Kondo_meas':
        arr = cu_extract(d, ['shifted_spectra'])
        for i, n in enumerate(arr):
            xi = np.array(n[0])
            m = da[i]['Kondo']['Kondo_window']
            mask = (xi >= m[0]) & (xi <= m[1])
            x.append(np.array(n[0])[mask])
            y.append(np.array(n[2])[mask])      
            
    elif spec == 'peak_0Vf' or spec == 'peak_0Vbf':
        arr = cu_extract(d, ['peak_0V','split_lorentz'])
        for i, n in enumerate(arr):
            x.append(np.array(n['V-ax']))
        if spec == 'peak_0Vf':
            for i, k in enumerate(arr):
                y.append(np.array(k['fit_1']))
        elif spec == 'peak_0Vbf':
            for i, k in enumerate(arr):
                y.append(np.array(k['best_fit']))

    return x, y

def propflt(dcti, key):
    tu = key_dct()
    if key in tu.keys():
        tu = tu[key]
        prop = tu[2]
        dct = dcti[tu[1][0]]   
        pr = []
        if len(prop) == 1:
            for i, k in enumerate(dct):
                pr.append(k[prop[0]])
        if len(prop) == 2:
            for i, k in enumerate(dct):
                pr.append(k[prop[0]][prop[1]])
        elif len(prop) == 3:
            for i, k in enumerate(dct):
                pr.append(k[prop[0]][prop[1]][prop[2]])        
        elif len(prop) == 4:
            for i, k in enumerate(dct):
                pr.append(k[prop[0]][prop[1]][prop[2]][prop[3]])  
    else:
        pr=dcti[key][0]
    return pr

def parachecker(di, para):
    nams = []
    for i,j in enumerate(di):
        if para not in j.keys():
            nams.append(j['name'])
    return nams

def prepare_curves(array, cut, norm):
    """
    norm is the used normalisation: 'ma' for max, 'mi' for min or None for
    no narmalisation
    """
    yuc = array[1]# y uncut
    if len(cut) > 0:
        V, y = cutter(array, cut)
        y_ = cutter(array, cut)[1] # need to make a second array indipendent from the first..
    else:
        V = array[0]
        y = array[1]
    Vlim, hiGn, loGn, vertsn, yn, ver = [], [], [], [], [], []
    if norm != None:
        yn_ = [] #yn.copy()
        for i, nam in enumerate(y):
            if norm == 'ma':
                yn.append(nam/max(nam))
                yn_.append(y_[i]/max(y_[i]))
            elif norm == 'mag': # mag --> maximum global
                yn.append(nam/max(nam))
                yn_.append(y_[i]/max(yuc[i]))
            elif norm == 'mi':
                yn.append(nam/min(nam))
                yn_.append(y_[i]/min(y_[i]))
            elif norm == 'mai':
                yn.append(nam/max(nam)-min(nam))
                yn_.append(y_[i]/max(y_[i]-min(y_[i])))
    elif norm == None:
        yn = y
        yn_ = y_
    #maxlen = max([len(i) for i in yn])
    
    carr = [V, yn_]
    for z in range(len(yn)):
        ysn = yn[z]
        xs = V[z]
        ys = yn_[z]
        ysn[0], ysn[-1] = 0, 0
        vertsn.append(list(zip(xs, ysn)))
        ver.append(list(zip(xs, ys)))
        hiGn.append(max(ysn))
        loGn.append(min(ysn))
        Vlim.append(max(abs(xs)))
        if z == 0:
            xes = V[z]
            ysta = yn_[z]
        else:
            xes = np.hstack([xes, V[z]])
            ysta = np.hstack([ysta, yn_[z]])
    extarr = [xes,ysta]            
    vertstn = np.transpose(vertsn)
    return [vertsn, vertstn, ver], hiGn, loGn, Vlim, carr 

def smooth(y,filterparams):#=[21,0,0,0.5,'interp']):
    """y can be a List of arrays"""
    window_sm = filterparams[0]
    polyorder = filterparams[1]
    deriv = filterparams[2]
    delta_sm = filterparams[3]
    mode = filterparams[4]
    Sm = []
    #for i in range(len(y)):
    Smi = [signal.savgol_filter(y, window_sm,\
                    polyorder, deriv, delta_sm, mode=mode)]
        #Sm.extend(Smi)
    return Smi

def smooth_param(x,sm_degree):
    """ degree of smoothin can be 'rough', 'middle', or 'fine'.
    """
    lx = len(x)
    if lx < 150:
        lx = lx + (150 - lx)
    if sm_degree == 'rough':
        win = round(lx/90)
        if win%2 == 0:
            win = win+1
        polyorder = 2
        deriv = 0
        delta = 0
        last = 'interp'
    elif sm_degree == 'middle':
        win = round(lx/40)
        if win%2 == 0:
            win= win+1
        polyorder = 4
        deriv = 0
        delta = 0
        last = 'interp'
    elif sm_degree == 'fine':
        win = round(lx/15)
        if win%2 == 0:
            win=win+1
        polyorder = 6
        deriv = 0
        delta = 0
        last = 'interp'
    return [win, polyorder, deriv, delta, last]

def cc(arlen):
    co = cols
    mcs = []
    for i in range(arlen):
        mc = mcolors.to_rgba(co[i], alpha=0.6)
        mcs.append(mc)
    return mcs

def cutter(arrs, cut):
    if len(cut) > 1:
        V, y = [], []
        #for i in range(len(arrs[0])):
        for n, j in enumerate(arrs[0]):
            ma = (j >= cut[0]) & (j <= cut[1])
            V.append(j[ma])
            y.append(arrs[1][n][ma])
        arr = V, y   # np.vstack([V,y,arrs[2]])
    else:
        arr = arrs
    return arr


def mkp(fp0,fp1): # makepath
    if os.path.exists(fp0) == True:
        pass
    elif os.path.exists(fp0) == False:
        os.makedirs(fp0)
        print('---path created--- \n'+fp0+' \n --- path created ---')
    fp = fp0+fp1
    return fp

def normto__(arr, norm='max'):
    if norm == 'max':
        marr = max(arr)
        narr = [i/marr for i in arr]
    elif norm == 'first':
        marr = arr[0]
        narr = [marr/i for i in arr]
    return narr

def new_dict1(filt_dicts, names, mk_):
    fd_ = filt_dicts
    na = {}
    na['swnr'] = (names, 'sweep_nr', 'sweep_number')
    #== new parameters ==========================================================
    #tupel On Off  
    oo = [propflt(fd_,'Gma'),propflt(fd_,'Gmi')]
    tuoo = ([i/j for i,j in zip(oo[0],oo[1])], '$G_{max} \ / \ G_{min}$')
    tuoo = tuoo + ('Conductance_On_Off_ratio',) # tuple on off
    na['G_onof'] = tuoo # new attribute dictionary
    #motorpos norm
    Mposn = normto__(propflt(fd_, 'Mpos'),'first')
    na['Mposn'] = (list(np.array(Mposn)*100), 'Mpos_rel [a.u.]', 'rel_Motorpos.')
    
    if mk_[0] == 'NDR':
        pAIndr=([abs(i) for i in propflt(fd_,'AIndr')],'$Area_{NDR} \ [A]$')
        na['pAIndr'] = pAIndr + ('NDR_Area_absolute',)
        
        pGIndr=([abs(i) for i in propflt(fd_,'GIndr')],#Integrated NDR area absolute value
                '$|G_{NDR\u222b}| \ [G/G_{0}]$')
        na['pGIndr'] = pGIndr + ('NDR_Integral_absolute',)
        
        ndr_r = [i-j for i,j in zip(propflt(fd_,'ndrcrr'), \
                                    propflt(fd_,'ndrcrl'))] # range (ndrcrr, ndrcrl are zero crossings)
        na['ndr_r'] = (ndr_r, 'NDR Range [mV]' , 'NDR Range')

    Imm = [i/abs(j) for i,j in zip(propflt(fd_,'Ima'),propflt(fd_,'Imi'))] # Imax/Imin
    na['Imm'] = (Imm, '$I_{max} \ / \ |I_{min}|$', 'Imax_by_Imin')
    
    return na

def new_dict2(dicts, pxy, rang, G_fac, names):
    filt  = filt_from_peaks(pxy, rang, G_fac) # returns: far, faral, farbi, index
    ndp = {}  # new attribute dictionary peaks
    ndp['Gpxr'] = (filt[0][0], '$U_{Bias} \ [mV]$', 'peak_position')
    ndp['Gpyr'] = (filt[0][1], '$ height \ [G/G_{0}]$', 'peak_height')
    ndp['swnr'] = (names, 'sweep_nr', 'sweep_number')
    Mposn = normto__(propflt(dicts, 'Mpos'),'first')
    ndp['Mposn'] = (list(np.array(Mposn)*10), 'Mpos_rel [a.u.]', \
                     'rel_Motorpos.')
    Gpn = [i/j for i,j in zip(filt[0][1],propflt(dicts,'Gma'))] #?? di y axes relative (normalized to max) ???
    ndp['Gpe_yn'] = (Gpn, '$rel \  height \ \  G_{peak} \  / \ |G_{max}|$ ',\
                  'rel_peak_height')
        
    return ndp

def plot_arrays(keys, dct, spdct):
    """

    keys : the keys to plot from standard or special dictionary
        DESCRIPTION.
    adct : dictionary
        standard dictionary (fids)
    spdct : dictionayry
        special dictionary

    Returns
    -------
    None

    """
    atdc =  key_dct() #attribute dictionary 
    vals, labels, ln = [], [], [] # values # ln : longname
    for i, j in enumerate(keys):
        if j in atdc.keys():
            vals.append(propflt(dct,j))
            labels.append(atdc[j][3])
            ln.append(atdc[j][0])
        else:
            vals.append(spdct[j][0])
            labels.append(spdct[j][1])
            ln.append(spdct[j][2])
    return [vals, labels, ln]


def plot_conc(xlst, ylst, zkey, dic, lonas):
    
    labels = [] # ln : longname
    tu = key_dct()
    if zkey != None:
        zar = propflt(dic,zkey) # z array
        labels = tu[xlst[0]][3], tu[ylst[0]][3], tu[zkey[0]][3]
    else:
        zar = np.zeros(len(dic[0]))
        labels = tu[xlst[0]][3], tu[ylst[0]][3]
    if len(xlst) > 1:
        zarr = zar + zar
        xst = np.hstack([propflt(dic,xlst[0]), \
                         propflt(dic,xlst[1])]).tolist() #x-stacked
        yst = np.hstack([propflt(dic,ylst[0]), \
                         propflt(dic,ylst[1])]).tolist() #y-stacked
        lnx = tu[xlst[0]][0]+'-'+tu[xlst[1]][0]
        lny = tu[ylst[0]][0]+'-'+tu[ylst[1]][0]
        if zkey != None:
            lnz  = tu[zkey][0]
    elif len(xlst) == 1:
        zarr = zar 
        xst = propflt(dic,xlst[0])
        yst = propflt(dic,ylst[0])
        lnx = tu[xlst[0]][0]
        lny = tu[ylst[0]][0]
        if zkey != None:
            lnz  = tu[zkey][0]
    x,y,z = [],[],[]
    for i,j in enumerate(xst):
        x = x+j
        y= y+yst[i]
        #if zkey != None:
        zlst = np.linspace(zarr[i],zarr[i],len(j)).tolist()
        z = z+zlst
    if zkey != None:
        vals =  [x,y,z]
        ln = [lnx,lny,lnz]
        if lonas != None:
            ln = lonas 
    else:
        vals = [x,y]
        ln = [lnx,lny]
    return [vals, labels, ln]
    
    
def filt_from_peaks(lstxy, lim, fac=None):
    lstx, lsty = lstxy
    farx, fary, index = [], [], [] #filtered array
    faral, farbi = [], []  # filtered array all, more than one entry
    for h, k in enumerate(lstx):
        k = np.array(k)
        if fac != None:
            y = np.array(lsty[h])/fac[h]
        else:
            y = np.array(lsty[h])
        flx = k[(k >= lim[0]) & (k <= lim[1])] #filtered list
        fly = y[(k >= lim[0]) & (k <= lim[1])]
        faral.append(flx)
        if len(flx) == 1:
            farx.append(flx[0])
            fary.append(fly[0])
            index.append(h)
        elif len(flx) > 1:
            farbi.append((flx,h))
    far = [farx,fary]
    return far, faral, farbi, index




def fits(x,y,model):
    params = model.guess(y, x)
    result = model.fit(y, params, x=x)
    comps = result.eval_components()
    report = result.fit_report(min_correl=0.5)
    best_fit = result.best_fit
    #gueses = (peak_positions, am)
    #print(result.fit_report(min_correl=0.5))
    return result, comps, best_fit, report, [], model, params

def Fit_linear(x, y):
    """other background could be quadratic or none"""
    #prefix = 'lin'
    shape= LinearModel()
    fit = fits(x, y, shape)
    return fit[0], fit[1], fit[2], fit[3], fit[4], fit[5], fit[6]


def savename(which,fig,what,val,comp,date):
    strval = str(val)
    if comp == '<':
        compstr = 'under'
    elif comp == '>':
        compstr = 'over'
    else:
        compstr = comp
    path = '/'+date + '/'+ which +'/'
    name = '/'+ what + compstr + strval+'_'+ fig +'.png'
    title = what+compstr + strval
    return name, title, path




# Figures 

def lineplot(lines, loG, hiG, Vlim, yla, cur, minz, negsw, azm = -60):

    edg = [] # list for edgecolors
    for i, j in enumerate(negsw):
        if j == 'n':
            edg.append('r')
        elif j == 'p':
            edg.append('k')
        elif j == 'z':
            edg.append('c')

    figl = plt.figure()
    ax = figl.add_subplot(1,1,1, projection='3d', azim=azm,proj_type= 'persp')
    #ax = fig.gca(projection='3d')

    #poly = PolyCollection(verts, edgecolors=[cc('r'), cc('g'), cc('b'), cc('y')])
    poly = LineCollection(lines, colors=cols, linewidths=(2,2,2,2))
    poly.set_alpha(0.7)

    ax.add_collection3d(poly, zs=cur, zdir='y')
    #PolyCollection
    ax.set_xlabel('Bias [mV]',fontsize = 14)
    ax.set_xlim3d(-max(Vlim),max(Vlim))
    ax.set_ylabel('sweep number',fontsize = 14)
    ax.set_ylim3d(min(cur), max(cur))
    ax.set_zlabel(yla,fontsize = 14)
    ax.set_zlim3d(minz, max(hiG))#
    ax.tick_params(axis='x' , labelsize = 12.0)
    ax.tick_params(axis='z' , labelsize = 12.0)
    ax.tick_params(axis='y' , labelsize = 11.0)
    figl.tight_layout()
    #plt.show()
    return figl

def edgecolor(param, ent):
    """
    ----------
    param : list
        parameters from measurement
    ent : tuple of 2 or 3 different entries
        entries in param (max 3 different)

    Returns: list of edgecolors
    -------
    None.

    """
    edg = [] # list for edgecolors
    for i, j in enumerate(param):
        if j == ent[1]:
            edg.append('r')
        elif j == ent[2]:
            edg.append('k')
        elif j == ent[3]:
            edg.append('c')
    return edg

def polyplot(Vets, loG, hiG, Vlim, yla, cur, minz, param, ent, azm = -60):
    """
    Parameters
    ----------
    Vets : TYPE
        DESCRIPTION.
    loG : TYPE
        DESCRIPTION.
    hiG : TYPE
        DESCRIPTION.
    Vlim : TYPE
        DESCRIPTION.
    yla : TYPE
        DESCRIPTION.
    cur : TYPE
        DESCRIPTION.
    minz : TYPE
        DESCRIPTION.
    param : TYPE
        DESCRIPTION.
    ent : TYPE
        DESCRIPTION.
    azm : TYPE, optional
        DESCRIPTION. The default is -60.

    Returns
    -------
    fig : TYPE
        DESCRIPTION.

    """
    edg = edgecolor(param, ent)
    
    plt.close('all')
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1, projection='3d',azim=azm,proj_type= 'persp')#'ortho')

    poly = PolyCollection(Vets,facecolors=cc(len(Vets)), edgecolors = edg, \
                          linewidths=2.0) #offsets = -1
    poly.set_alpha(0.4)
    ax.add_collection3d(poly, zs=cur, zdir='y')
    #PolyCollection
    ax.set_xlabel('Bias [mV]',fontsize = 14)
    ax.set_xlim3d(-max(Vlim),max(Vlim))
    ax.set_ylabel('sweep number',fontsize = 14)
    ax.set_ylim3d(min(cur), max(cur))
    ax.set_zlabel(yla,fontsize = 14)
    ax.set_zlim3d(minz, max(hiG))#
    ax.tick_params(axis='both' , labelsize = 12.0)
    fig.tight_layout()
    plt.show()
    #plt.close(fig)
    return fig

def fig1_multi(V,y, ylab, legs, xlab='Bias [mV]', bf = 12,axvl=[]):
    """
    signal 'iets' or 'dGdV', I or didv
    """
    co = cols
    #com = cm.tab20c
    #plt.gca().set_prop_cycle(plt.cycler('color', com(np.linspace(0, 1, len(y)))))
    #lab = np.arange(0,10,1)
    #I,LIA1,LIA2, Iten = y[0], y[1], y[2], y[3]
    f1 = plt.figure(3,figsize=[7.2, 4.8])
    ax1 = f1.add_subplot(111)
    for i, j in enumerate(y):
        #co_ = colormap(i) 
        ax1.plot(V[i], j, co[i], label='sw ' + str(legs[i]) )
        #ax1.plot(V[i], j, label='sw ' + str(legs[i]) )
    #ax4 = ax3.twinx()
    f1.tight_layout()
    ax1.set_ylabel(ylab, fontsize = 16)
    #ax1.plot(V, LIA1, c='K')
    if len(axvl)>0:
        for i,j in enumerate(axvl):
            colo = ['r', 'r', 'k', 'k', 'b', 'b']
            ax1.axvline(j, linestyle = '--', c = colo[i], linewidth= 0.9)
    ax1.set_xlabel(xlab, fontsize = 16)
    ax1.tick_params(axis='both' , labelsize = 14.0)
    params = {'legend.markerscale':15}
    plt.rcParams.update(params)
    ax1.legend(fontsize=12)#,loc = 'lower left')
    ax1.legend(bbox_to_anchor=(1,1), fontsize = bf, loc = 'upper left')
    f1.tight_layout()
    return f1

def fig_stack (arr,ylab,legs,ti,xlab='Bias [mV]',ydist=0.1,hgt=7.8,bf=12):
    """
    signal 'iets' or 'dGdV', I or didv
    ydist --> distance between plots
    hgt: height
    """
    V, y = arr
    co = cols
    #lab = np.arange(0,10,1)
    #I,LIA1,LIA2, Iten = y[0], y[1], y[2], y[3]
    f1 = plt.figure(3,figsize=[4.8, hgt])
    ax1 = f1.add_subplot(111)
    for i in range(len(y)):
        ax1.plot(V[i], y[i]+i*ydist, co[i], label='sw ' + str(legs[i]))
        y0 = np.zeros((len(V[i]), 1))+i*ydist
        ax1.plot(V[i], y0 , '--', color = co[i], linewidth=0.8 )
    #ax4 = ax3.twinx()
    f1.tight_layout()
    ax1.set_ylabel(ylab, fontsize = 16)
    ax1.set_title(ti, size=16)
    #ax1.plot(V, LIA1, c='K')
    ax1.set_xlabel(xlab, fontsize = 16)
    plt.grid(True, which='major',axis='x',color='g',linewidth=0.4)
    ax1.axes.yaxis.set_ticks([])
    ax1.tick_params(axis='x', labelsize = 14.0)
    params = {'legend.markerscale':15}
    plt.rcParams.update(params)
    ax1.legend(bbox_to_anchor=(1,1), fontsize = bf, loc = 'upper left')
    f1.tight_layout()
    return f1

colors = [mcolors.to_rgba(c)
          for c in plt.rcParams['axes.prop_cycle'].by_key()['color']]
    

def figmulti_sc(arr, si, tis, xln,yln,da,pa, sav = False, dpi = 250):
    """ signal 'iets' or 'dG/dV' """
    dire, mfld = pa
    x1,y1,x2,y2,x3,y3 = arr[0]
    xl1,yl1,xl2,yl2,xl3,yl3 = arr[1]
    lns = arr[2] #long names
    if len(tis) == 0:
        tis = np.array(lns)#[1::2]
    plt.close('all')
    f1 = plt.figure(1,figsize=(6.4,7.6))
    ax1 = f1.add_subplot(212)#, sharex=ax3)
    ax2 = f1.add_subplot(222)
    ax3 = f1.add_subplot(221) #, sharex=ax3)
    #ax3.margins(x=0, y=-0.25)
# =============================================================================
#     ax3 = f1.add_subplot(313)
#     ax1 = f1.add_subplot(311, sharex=ax3)
#     ax2 = f1.add_subplot(312, sharex=ax3)
# =============================================================================
    #ax4 = ax3.twinx()
    f1.tight_layout()
    ax1.scatter(x1,y1, c='k')
    ax1.set_ylabel(yl1, fontsize = 16)
    ax1.set_xlabel(xl1, fontsize = 16)
    ax1.set_title(tis[1], fontsize = 12)
    ax1.tick_params(axis='both' , labelsize = 12.0)
    
    ax2.scatter(x2,y2, c='k', s=si)
    ax2.set_ylabel(yl2, fontsize = 16)
    ax2.set_xlabel(xl1, fontsize = 16)
    ax2.set_title(tis[3], fontsize = 12)
    ax2.tick_params(axis='both' , labelsize = 12.0)
    
    #lower big plot
    ax3.scatter(x3, y3, c='k', s=si)#, label='dG/dV [G0]')
    ax3.set_ylabel(yl3, fontsize = 16)
    ax3.set_xlabel(xl3, fontsize = 16)
    ax3.set_title(tis[5], fontsize = 12)
    ax3.tick_params(axis='both' , labelsize = 12.0)
    #ax3.legend(['T = ' + str(temp)+ ' K'])
    f1.tight_layout()
    if sav == True :
        svn = savename(mfld,'multi',xln,yln,'_vs_',da)
        sap = mkp(dire+ '/Evaluation_plots/'+svn[2], svn[0])
        f1.savefig(sap,dpi=dpi,format='jpg')
    return f1   

def f_scat(vals, ti,labs,design='bone',si=30, xlock ='yes',ano=[],fit='n'):
    """ 
    yl : y label
    xl : x label
    bl : bar label
    ano: annotate
    """
    V,yv,param = vals
    xl,yl,bl = labs
    plt.close('all')
    f1 = plt.figure(3,figsize=[6.4, 4.8])
    ax1 = f1.add_subplot(111)
    
    if len(ano) > 0:
        txt = [str(i) for i in ano]
        for lab, x, y in zip(txt, V, yv):
            ax1.annotate(lab, xy = (x, y), xytext=(x,y), color = 'darkgray', \
                      ha='center',va='center', \
                          fontsize = 6)
    if fit == 'y':
        lf = Fit_linear(V,yv)[0].best_values
        lin = np.array(V)*lf['slope'] + lf['intercept']
        lin=lin.tolist()
        lfr= str(round(lf['slope'],2)),str(round(lf['intercept'],2)) # IV fit results
        ti = 'y = ' + lfr[0] +' * x' + ' +- '+ lfr[1]
        ax1.plot(V, lin, ':', c='r', lw=1.2)
        
    plot = ax1.scatter(V, yv, c=param, cmap=design, s=si, alpha=0.75)
    cbar = f1.colorbar(plot)
    cbar.set_label(bl, size = 14)
    #f1.tight_layout()
    ax1.set_title(ti,size=16)
    if xlock == 'yes':
        plt.xticks(V[::2])
        ax1.grid(which='major',axis='x',color='k',linestyle=':',\
                 linewidth=0.5, alpha=0.85)
    #f1.subtitle('OV Bias Conductivity < 0.05 G0',fontsize=16)
    ax1.set_ylabel(yl, fontsize = 16)
    ax1.set_xlabel(xl, fontsize = 16)
    ax1.tick_params(axis='both',labelsize=14)#,labelbottom='off',labeltop='on')
    ax1.legend([plot],['N = ' + str(len(yv))+ ' points'], markerscale=1)
    #plt.legend([plot], [legend],  fontsize = 12, loc= lo, markerscale=10)
    f1.tight_layout()
    return f1

def f_scat_nb(vals,ti,labs,design='bone',si=32,xlock='yes',ano=[],fit='n'): # figure scatter no bar (colorbar)
    """ 
    yl : y label
    xl : x label
    bl : bar label
    """
    V,yv = vals
    xl, yl = labs
    plt.close('all')
    f1 = plt.figure(3,figsize=[5.6, 4.8])
    ax1 = f1.add_subplot(111)
    plot = ax1.scatter(V, yv, c='k', s=si, alpha=0.85)
    if len(ano) > 0:
        txt = [str(i) for i in ano]
        for lab, x, y in zip(txt, V, yv):
            ax1.annotate(lab, xy = (x, y), xytext=(x,y), color = 'silver', \
                      ha='center',va='center', \
                          fontsize = 6)
    if fit == 'y':
        lf = Fit_linear(V,yv)[0].best_values
        lin = np.array(V)*lf['slope'] + lf['intercept']
        lin=lin.tolist()
        lfr= str(round(lf['slope'],2)),str(round(lf['intercept'],2)) # IV fit results
        ti = ['Y = ' + lfr[0] +' *X' + ' +- '+ lfr[1]]
        ax1.plot(V, lin, ':', c='r', lw=1)
    if xlock == 'yes':
        plt.xticks(V[::2])
    ax1.set_title(ti,size=16)
    ax1.set_ylabel(yl, fontsize = 16)
    ax1.set_xlabel(xl, fontsize = 16)
    ax1.tick_params(axis='both',labelsize=14)#,labelbottom='off',labeltop='on')
    ax1.grid(which='major',axis='x',color='k',linestyle=':',\
             linewidth=0.5, alpha=0.85)

    ax1.legend(['N = ' + str(len(yv))+ ' points'], markerscale=1)
    f1.tight_layout()
    return f1

def fig_ps(arr,mfld,dire,da,ti =' ',sty='jet',sav=False,si=45,dpi=200,\
           ano=[], fit = 'n', hist='no', logy= 'no', bins='50'):         # figure properties scatter
    """
    Parameters
    ----------
    arr : list
        3 lists [values, labels, longnames] for x- y- axis and colorbar (if)
    mfld : string
        mainfolder
    dire: string
          Directory where fig should be saved
    da : string
        date 
    cbar : string, optional
        plot with colorbar or not. The default is 'yes'.
    ti : string, optional
        Title of the plot. The default is [].
    sty : string, optional
        colorbar colormap. The default is 'jet'.
    sav : Boolean 
        save the plot --> True or False
    si : int
        size of the scatters
    dpi: size of the saved figure in dpi
    Returns 
    Figure and save that (if sav == 'yes')
    """
    #xv = propflt(dct[x[1]], x[2])
    if len(arr[0]) == 3:
        xln, yln,bln = arr[2] # ln : long name
        if arr[1][0] == 'sweep_nr':
            lock = 'yes'
        else: 
            lock = 'no'
        svn = savename(mfld,bln+'_'+sty,xln,yln,'_vs_',da)
        if hist=='no':
            f = f_scat(arr[0], ti, arr[1],sty, si, lock=lock,ano=ano,fit=fit) 
        else:
            f = fig_2dHist(arr[0], ti, arr[1], bins, sty, logy, trans = 'no')
    elif len(arr[0]) == 2:
        x,y = arr[0]
        xl,yl = arr[1]
        xln, yln = arr[2] # ln : long name
        if xl == 'sweep_nr':
            lock = 'yes'
        else: 
            lock = 'no'
        svn = savename(mfld,sty,xln,yln,'_vs_',da) 

        f = f_scat_nb(arr[0], ti, arr[1],'bw', si, xlock=lock,ano=ano,fit=fit)   
    sap = mkp(dire+ '/Evaluation_plots/'+svn[2], svn[0])
    if sav == True :
        f.savefig(sap,dpi=dpi,format='jpg')
    return f


def fig_2dHist(data, Title,lab, bins, comap, logy, trans = 'yes'):
    if trans == 'yes':
        x,y = arr_trans(data)
    else:
        x,y = data
    xl,yl = lab
    xmin = x.min() #np.log10(x.min())
    xmax = x.max() #np.log10(x.max())
    if logy == 'yes':
        ymin = np.log10(y.min())# - np.log10(y.min())/10
        ymax = np.log10(y.max())# + np.log10(y.max())/10
        ybins = np.logspace(ymin, ymax, bins) # <- make a range from 10**ymin to 10**ymax
        
    elif logy == 'no':
        ymin = y.min()# - np.log10(y.min())/10
        ymax = y.max()# + np.log10(y.max())/10
        ybins = np.linspace(ymin, ymax, bins) # <- make a range from 10**ymin to 10**ymax
    
    xbins = np.linspace(xmin, xmax, bins) # <- make a range from 10**xmin to 10**xmax

    hdat = np.histogram2d(x, y, bins=(xbins, ybins))
    counts, _, _ = np.histogram2d(x, y, bins=(xbins, ybins))
    
    fig = plt.figure(3,figsize=(11, 5.8))
    axes  = fig.add_subplot(111)
    pcm = axes.pcolormesh(xbins, ybins, counts.T, cmap=comap)
    #pcm = axes.pcolormesh(hdat[1], hdat[2], counts.T, cmap=comap)
    plt.colorbar(pcm)
    #fig, axes = plt.subplots(nrows=2, ncols=2)
    axes.set_title(Title, fontsize = 24)  
    #axes.hist2d(data[:, 0], data[:, 1], bins=bins)#, norm=mcolors.LogNorm())#, label)
    
    #plt.ylim(bottom = min(data[:, 1])-min(data[:, 1])/2)
    #plt.ylim(top=max(data[:, 1])+ max(data[:, 1])/2)
    if logy == 'yes':
        axes.set_yscale('log')
    axes.set_ylabel(yl, fontsize = 20)
    axes.set_xlabel(xl, fontsize = 20)
    axes.tick_params(axis='both' , labelsize = 18.0)
    
    axes.set_xlim(xmin=xbins[0])
    axes.set_xlim(xmax=xbins[-1])
    axes.set_ylim(ymin=ybins[0])
    axes.set_ylim(ymax=ybins[-1])

    fig.tight_layout()
    return fig

def arr_trans(arr):
    for i,j in enumerate(arr[0]):
        x = j
        y = arr[1][i]
        if i == 0:
            ext_x = j # extended x axis
            ext_y = y
        else:
            ext_x = np.hstack((ext_x, x))
            ext_y = np.hstack((ext_y, y))
    return ext_x, ext_y

def fig_average(arr, bins, xlab, ylab, logy='no', sm='middle', axvl=[]):
    x,y = arr_trans(arr)
    xmin = x.min() 
    xmax = x.max() 
    if logy == 'yes':
        ymin = np.log10(y.min())
        ymax = np.log10(y.max())
        ybins = np.logspace(ymin, ymax, bins) # <- make a range from 10**ymin to 10**ymax
        
    elif logy == 'no':
        ymin = y.min()
        ymax = y.max()
        ybins = np.linspace(ymin, ymax, bins) # <- make a range from 10**ymin to 10**ymax
    xbins = np.linspace(xmin, xmax, bins) # <- make a range from 10**xmin to 10**xmax
    hdat = np.histogram2d(x, y, bins=(xbins, ybins))
    av_, sums = [], []
    for i,j in enumerate(hdat[0]):
        Gm = hdat[2][:-1]*np.transpose(j)
        su = np.sum(j)
        if i == 0:
            Gmatrix = Gm
        else:
            Gmatrix = np.vstack([Gmatrix, Gm])
        sums.append(su)

    hdt = hdat[1][:-1]
    for i,j in enumerate(Gmatrix):
        ave = np.sum(j)/sums[i]
        av_.append(ave)
    if sm == 'rough' or sm == 'middle' or sm == 'fine':
        av = smooth(av_, smooth_param(av_, sm))[0]
    elif sm == '' or sm == 'no':
        av = av_
    # make figure:
    plt.close('all')
    f1 = plt.figure(3,figsize=[7.2, 4.8])
    ax1 = f1.add_subplot(111) 
    ax1.plot(hdt, av, color = 'k')
    f1.tight_layout()
    ax1.set_ylabel(ylab, fontsize = 16)
    if len(axvl)>0:
        for i,j in enumerate(axvl):
            colo = ['r', 'r', 'k', 'k', 'b', 'b']
            ax1.axvline(j, linestyle = '--', c = colo[i], linewidth= 0.9)
    ax1.set_xlabel(xlab, fontsize = 16)
# =============================================================================
#     ax1.tick_params(axis='both' , labelsize = 14.0)
#     params = {'legend.markerscale':15}
#     plt.rcParams.update(params)
#     ax1.legend(fontsize=12)#,loc = 'lower left')
#     ax1.legend(bbox_to_anchor=(1,1), fontsize = bf, loc = 'upper left')
# =============================================================================
    f1.tight_layout()
    return f1





############################ explainer ####################################33

def explainer(section):
    expl = {}
    expl['0'] = 'load modules needed in the programm. The fcts_summarizer\
                      \n exhibits functions for raw data handling and plotting \
                      \n written by ourselfs'
    expl['1'] = 'define the directory where the folders curve_attributes, \
                    \n curves_analysed, and extrema_analysed are. The console \
                        \n tells you if that directory is existing.'
    expl['2'] = 'define the number of curves you want to consider. \
                      \n that is either done by writing every number in a \
                          \n dictionary or by writing the first and last one \
                              and let numpy create an array. '
                        
    
    return print(expl[section])
