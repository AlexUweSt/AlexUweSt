#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 06:48:07 2021

@author: filipk
functions for the ACP_XXX.py programm 
ACP --> auto curve plotter

"""

import json
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
#from numpy import diff
#import pandas as pd
import os
#from pathlib import Path
import matplotlib as mpl
import math
import time
import scipy.constants as spc
from scipy import signal
import scipy.integrate as integ
import datetime as dt
import shutil as sht
import copy
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)


def filefinder(directory):
    a = []
    for file in os.listdir(directory):
        a.append(file)
    b = [k for k in a if "IETS" in k]  # or "upsweep"or "downsweep"
    c = [k for k in a if "upsweep" in k]
    d = [k for k in a if "downsweep" in k]
    e = [k for k in a if "G0" in k]
    j = [str(directory+x) for x in b]   # fullsweep Matrix paths
    # firstsweep ('upsweep' - Matrix) paths
    l = [str(directory+x) for x in c]
    # secondsweep ('downsweep' - Matrix) paths
    m = [str(directory+x) for x in d]
    o = [str(directory+x) for x in e]   # GO paths
    print('    sweeps: ' + str(len(j)) + '    OCs: ' + str(len(o)))
    return [b, c, d, e], [j, l, m, o]


def filesorter(fp):
    fpf = fp[0]  # filepath fullsweep
    fp1 = fp[1]
    fp2 = fp[2]
    fulls, firs, secs, ocs = [], [], [], []  # fulsort, firstsort, ..
    fullsp, firsp, secsp, ocsp = [], [], [], []
    funams, finams, senams, ocnams = [], [], [], []
    for i, j in enumerate(fpf):
        funam = os.path.split(j)[1]
        funams.append(funam)
        indf = funam.find('IETS0')  # fullsweepname
        fulls.append(int(funam[indf+4:indf+8]))
        fullsp.append(j)

        finam = os.path.split(fp1[i])[1]
        finams.append(finam)
        indfi = finam.find('upsweep0')  # fullsweepname
        if indfi >= 0 :  
            firs.append(int(finam[indfi+7:indfi+11]))
        else:
            indfi = finam.find('downsweep0')
            firs.append(int(finam[indfi+9:indfi+13]))
        firsp.append(fp1[i])

        senam = os.path.split(fp2[i])[1]
        senams.append(senam)
        inds = senam.find('downsweep0')  # fullsweepname
        if inds >= 0:
            secs.append(int(senam[inds+9:inds+13]))
        else:
            inds = senam.find('upsweep0')
            secs.append(int(senam[inds+7:inds+11]))
        secsp.append(fp2[i])

    for i, j in enumerate(fp[3]):
        ocnam = os.path.split(j)[1]
        ocnams.append(ocnam)
        indo = ocnam.find('G')  # fullsweepname
        ocs.append(int(ocnam[indo+1:indo+5]))
        ocsp.append(j)

    xa = np.argsort(fulls)
    fuu = np.array(fullsp)[xa]
    funams_ = np.array(funams)[xa]

    xb = np.argsort(firs)
    fii = np.array(firsp)[xb]
    finams_ = np.array(finams)[xb]

    xc = np.argsort(secs)
    see = np.array(secsp)[xc]
    senams_ = np.array(senams)[xc]

    xd = np.argsort(ocs)
    occ = np.array(ocsp)[xd]
    ocnams_ = np.array(ocnams)[xd]

    return [fuu, fii, see, occ], [funams_, finams_, senams_, ocnams_]

def load_mainsweep(fipa, ad_info, ind):
    """
    filepath skiprows measurement program state (version)
    """
    state = ad_info['state']
    Meas_d = {}
    fisw_d = {}
    secsw_d = {}
    sr = ad_info['skiprows']
    if state == 5 or state == 6:
        MEAS = np.genfromtxt(fipa[0][ind], skip_header=sr, delimiter='\t', \
                         usecols=(0, 1, 2, 3, 4, 5, 6, 7, 8), dtype=np.str)
        times = np.genfromtxt(fipa[0][ind], skip_header=sr, delimiter='\t',
                              usecols=9, dtype='str')
        times = np.char.replace(times, '.', ' ')
        times = np.char.replace(times, ',', '.')
        times = [datetime.strptime(i, ' %H:%M:%S.%f %d %m %Y') for i in times]
    else:
        MEAS = np.loadtxt(fipa[0][ind], skiprows=sr, delimiter='\t', \
                          dtype=np.str)
    MEAS = np.char.replace(MEAS, ',', '.')
    MEAS = np.char.replace(MEAS, '\'', '')
    MEAS = np.char.replace(MEAS, 'b', '').astype(np.float64)
    MEAS_T = np.transpose(MEAS)
    
    if state == 6:
        uc = (0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
        header = np.genfromtxt(fipa[0][ind], dtype='str',skip_header=sr-1,\
                               max_rows=sr-(sr-1))#,usecols=uc)
        header = np.concatenate([header, ['time']])
        header_ = np.genfromtxt(fipa[0][ind], dtype ='str', max_rows=sr-2)
        header_ = np.char.replace(header_, ',', '.')
        header = np.vstack([header_, header])
    elif state == 5:
        uc = (0, 1, 2, 3, 4, 5, 6, 7, 8)
        header = np.genfromtxt(fipa[0][ind], dtype='str',skip_header=sr-3,\
                               max_rows=3,usecols=uc)
        #header = np.concatenate([header, ['time']])
        header_ = np.genfromtxt(fipa[0][ind], dtype ='str', max_rows=2, \
                                usecols=uc)
        header_ = np.char.replace(header_, ',', '.')
        header = np.vstack([header_, header])
    else:
        header = np.genfromtxt(fipa[0][ind], dtype ='str', max_rows=5)
        header = np.char.replace(header, ',', '.')
        
    if os.path.getmtime(fipa[1][ind]) < os.path.getmtime(fipa[2][ind]):
        fisw = load_updown(fipa[2][ind])
        secsw = load_updown(fipa[1][ind])
    else:
        fisw = load_updown(fipa[1][ind])
        secsw = load_updown(fipa[2][ind])

    if state == 6:
        for i, j in enumerate(header[4][:-1]):
            Meas_d[j] = MEAS_T[i]
            fisw_d[j] = fisw[i]
            secsw_d[j] = secsw[i]
        # = header[4][:-1]
    else:
        for i, j in enumerate(header[4]):
            Meas_d[j] = MEAS_T[i]
            fisw_d[j] = fisw[i]
            secsw_d[j] = secsw[i]
    if 'RL_LIA2' in Meas_d.keys():
        Meas_d['R_LIA2'] = Meas_d['RL_LIA2']
        fisw_d['R_LIA2'] = fisw_d['RL_LIA2']
        secsw_d['R_LIA2'] = secsw_d['RL_LIA2']
        header[4][2] = 'R_LIA2'
        del Meas_d['RL_LIA2']
        del fisw_d['RL_LIA2']
        del secsw_d['RL_LIA2']
    if state == 1:
        Meas_d['Bias_DMM'] = Meas_d['LIA_1_DMM']
        fisw_d['Bias_DMM'] = fisw_d['LIA_1_DMM']
        secsw_d['Bias_DMM'] = secsw_d['LIA_1_DMM]
        header[4][6] = 'Bias_DMM'
        del Meas_d['LIA_1_DMM']
        del fisw['LIA_1_DMM']
        del secsw['LIA_1_DMM']
        mh = header[4]
    elif state == 2:
        mh = header[4]
        Meas_d['Bias_DMM'] = MEAS_T[0]*2.0455-0.0005
        fisw_d['Bias_DMM'] = fisw[0]*2.0455-0.0005
        secsw_d['Bias_DMM'] = secsw[0]*2.0455-0.0005
        mh = np.append(mh, 'Bias_DMM')
    elif state == 3:
        mh = header[4]
        Meas_d['Bias_DMM'] = MEAS_T[0]*2.0455-0.0005
        fisw_d['Bias_DMM'] = fisw[0]*2.0455-0.0005
        secsw_d['Bias_DMM'] = secsw[0]*2.0455-0.0005
        mh = np.append(mh, 'Bias_DMM')
# =============================================================================
#     elif state == 4 or state == 5 or state == 6:
#         mh = header[4]
#         mh.append('time')
# =============================================================================
    mhead = []  
    for i, j in enumerate(Meas_d.keys()):
        mhead.append(j)
        if i == 0:
            meas_new = Meas_d[j]
            finew = fisw_d[j]
            secnew = secsw_d[j]
        else:
            meas_new = np.vstack((meas_new, Meas_d[j]))
            finew = np.vstack((finew, fisw_d[j]))
            secnew = np.vstack((secnew, secsw_d[j]))
            
    if state == 5 or state == 6:
        Meas_d['time'] = times
        np.vstack([meas_new,times])
        mhead.append('time')
                
    Meas_save = np.transpose(meas_new)
    fisave = np.transpose(finew)
    secsave = np.transpose(secnew)
    arry = [Meas_d, Meas_save, mhead]
    return arry, [fisw_d, fisave], [secsw_d, secsave], header


def load_updown(fipa):
    if os.path.getsize(fipa) < 1500:
        upsw = np.transpose(np.vstack((np.zeros(9), np.zeros(9))))
    else:
        upsweep = np.loadtxt(fipa, skiprows=0, delimiter='\t', dtype=np.str)
        upsweep = np.char.replace(upsweep, ',', '.')
        upsweep = np.char.replace(upsweep, '\'', '')
        upsweep = np.char.replace(upsweep, 'b', '').astype(np.float64)
        upsw = np.transpose(upsweep)
    return upsw


def Header(header, adi):
    """
    Parameters
    ----------
    Header : Header array
    state : number int or float
        MEasurement programm differentiation 
    adi : list
        aditional information that is not given in the header but known
    Returns
    -------
    new header and curvve infos needed for further treading the curves

    """
    Hdi = {}
    ty  = [int, int, int, int, float, int, float, \
                            int, int,  str]
    h3 = []
    for i,j in enumerate(header[3]):
        h3.append(j.astype(ty[i]))
        #h3.append(ty[i](j))
    #param = np.vstack(([header[1].astype(np.float)], [np.array(h3)]))
    if header[2][6] != 'LockAmpl':
        header[2][6] = 'LockAmpl'

    for i, j in enumerate(header[0]):
        Hdi[j] = float(header[1][i])
        Hdi[header[2][i]] = h3[i].item()

    if Hdi['LockAmpl'] == 0:
        lockampl = adi['LockAmpl']
        header[3][6] = lockampl
    elif 1 < Hdi['LockAmpl'] < 60:
        lockampl = Hdi['LockAmpl']/1E3
        #Hdi['LockAmpl'] = lockampl
        header[3][6] = lockampl
    elif 1E-6 < Hdi['LockAmpl'] < 0.9:
        pass
    else:
        lockampl = Hdi['LockAmpl']/1E5
        header[3][6] = lockampl
    Hdi['LockAmpl'] = lockampl 

    return Hdi, header


def load_oc(path, adic):
    state = adic['state']
    OC_file = np.loadtxt(path, delimiter='\t', dtype=np.str)
    OC_file = np.char.replace(OC_file, ',', '.')
    OC_file = np.char.replace(OC_file, 'b', '').astype(np.float64)
    oc = np.transpose(OC_file)
    arr = {}
    timestamp = arr['timestamp'] = oc[0]
    MposOC = arr['MposOC'] = oc[1]/1E+7
    vel = arr['vel'] = oc[2]
    Gval_ = arr['Gval_'] = oc[3]
    Rval = arr['Rval'] = oc[4]
    if state == 4 or state == 5:
        LIAG_ = arr['LIAG_'] = oc[7]
        #tempe = oc[11]
        LIA1Dmm = arr['LIA1Dmm'] = oc[5]
    if state == 6:
        LIAG_ = arr['LIAG_'] = oc[7]
        tempe = arr['tempe'] = oc[11]
        LIA1Dmm = arr['LIA1Dmm'] = oc[5]
    #LIA2Dmm = oc[6]
    return arr

def calc_voc(ocfile):
    G0 = spc.physical_constants['conductance quantum'][0]
    voc = 1/((1/(G*G0)+Rpre))*Vin*CAg
    return voc

def negvals(dct):
    Gval_ = dct['Gval_']
    LIAG_ = dct['LIAG_']
    Gval, LIAG = [], []
    for i, g in enumerate(Gval_):
        if g < 0:
            Gval.append(((g*(-1)+3)))
        else:
            Gval.append((g+0))
    for i, g in enumerate(LIAG_):
        if g < 0:
            LIAG.append(((g*(-1)+3)))
        else:
            LIAG.append((g))
    dct['LIAG_'] = LIAG
    dct['Gval'] = Gval

    return dct


def theta_(X, R, thetsh):
    """
    other method is theta = ((R**2 - X**2)**0.5)/X - Theta_ref 
    Parameters
    ----------
    X : list
        LIA x measurement
    R : list
        LIA R value
    thetsh : int or float
        thehta shift

    Returns
    -------
    thet : np array
        thetha
    thetm : np array
        thetha shiftet

    """
    thet = []
    tt = []
    for i, j in enumerate(R):
        if j <= 0:
            t = 1E-7
        else:
            t = j
        tt.append(X[i]/t)
    #tt= X/R
    for i in range(len(tt)):
        if tt[i] > 1:
            ttt = 1
        elif tt[i] < -1:
            ttt = -1
        else:
            ttt = tt[i]
        thet.append(math.acos(ttt))
    thet = np.array(thet)
    thetm = thet - thetsh
    return [thetm, thet]


def calc_arrays(Meas_dict, Head_dict, adi):
    """
    Parameters
    ----------
    Measarr : Dictionary
        Masurement array 
    Head_dict : dictionary
        Header Information in Dictionary form
    adi : dictionary
        aditional Information

    Returns
    -------
    Dictionary of transformed arrays of the measurements file

    """
    G0 = spc.physical_constants['conductance quantum'][0]
    md = Meas_dict
    hdic = Head_dict  # hdic --> Header dictionary
    R_DUT = adi['Rpre']
    R_pre_G = (1/R_DUT)/G0
    arrdct = {}
    Lampl = hdic['LockAmpl']
    CAg = hdic['CA_gain']
    #state = adi['state']

    LIA1_gain = hdic['LIASens1']/(hdic['LIAgain1']*10) * adi['facs'][1]
    LIA2_gain = hdic['LIAsens2']/(hdic['LIAgain2']*10) * adi['facs'][2]

    # Bias:
    arrdct['Bias_Vspt'] = md['Bias[V]']
    V = arrdct['Bias_DMM'] = md['Bias_DMM'] * 1000
    # voltagedrop over preresistor and sample bias
    us = arrdct['Usample'] = (V-(md['I(DMM)']*CAg*R_DUT))
    # Current
    arrdct['I(DMM)'] = us/(us/(md['I(DMM)'] * CAg)) * (1E6)*adi['facs'][0]
    # LIA_direct
    R_LIA_1a = arrdct['R_LIA1'] = md['R_LIA1']*CAg/(Lampl*G0)
    R_LIA_1 = arrdct['R_LIA_1corr'] = 1/((1/R_LIA_1a) - (1/R_pre_G))
    
    X_LIA2 = adi['facs'][2] * md['X_LIA2']*CAg/(G0*Lampl**2)
    arrdct['X_LIA2'] = X_LIA2
    
    if 'Th_LIA2' in md.keys():
        if 'lockin_phase' in hdic.keys():
            arrdct['Th_LIA_2'] = md['Th_LIA2'] - hdic['lockin_phase']
        else:
        # Theta LIA 2 SR830 outputs 1V per 18 degrees
            arrdct['Th_LIA_2'] = md['Th_LIA2']# * 18

    else:
        R_LIA2 = arrdct['R_LIA2'] = md['R_LIA2']*CAg/(G0*Lampl**2)
        if 'lockin_phase' in hdic.keys():
            arrdct['Th_LIA2']=theta_(X_LIA2,R_LIA2,hdic['lockin_phase'])[0]
        else:
            arrdct['Th_LIA2'] = theta_(X_LIA2, R_LIA2, 0)[0]
            
    X_LIA_1a = adi['facs'][1]*md['X_LIA1']*CAg/(Lampl*G0)
    arrdct['X_LIA1'] = X_LIA_1a
    X_LIA_1 = arrdct['X_LIA1corr'] = 1/((1/X_LIA_1a) - (1/R_pre_G))

    LIA2 = arrdct['LIA2_DMM'] = md['LIA_2_DMM']*CAg/(G0*Lampl**2)*LIA2_gain
    if 'LIA_1_DMM' in md.keys():
        LIA1_a = arrdct['LIA1_DMM'] = md['LIA_1_DMM']*CAg/(Lampl*G0)*LIA1_gain
        LIA1 = arrdct['LIA1corr'] = 1/((1/LIA1_a) - (1/R_pre_G))
        arrdct['XY'] = LIA2/LIA1
        res_Lia1 = arrdct['res_Lia1'] = (1/(LIA1*G0))/1000

    else:
        arrdct['XY'] = LIA2/X_LIA_1
        res_Lia1 = (1/(X_LIA_1*G0))/1000

    if 'Lia1_phase' in hdic.keys():
        arrdct['thet_LIA1'] = theta_(X_LIA_1, R_LIA_1, hdic['Lia1_phase'])[0]
    else:
        arrdct['thet_LIA1'] = theta_(X_LIA_1, R_LIA_1, 0)[0]


    return arrdct


def G_zero(meas_dict, hdic):
    """
    Parameters
    ----------
    meas_dict : dictionary
        measurement dictionary
    hdic : dictionary
        header dictionary

    Returns
    -------
    header dict with G_zero

    """
    if 'LIA_1' in meas_dict.keys() or 'LIA_1_DMM' in meas_dict.keys():
        z_len = int(np.round(len(meas_dict['LIA1corr'])/2))
        hdic['G_zero'] = meas_dict['LIA1corr'][z_len]
    else:
        z_len = int(np.round(len(meas_dict['X_LIA1corr'])/2))
        hdic['G_zero'] = meas_dict['X_LIA1corr'][z_len]
    return hdic


def G_scale(Dictionaries, hedict):
    Fuldi, fi_di, sec_di = Dictionaries
    G_keys = ['R_LIA_1', 'R_LIA_1corr', 'X_LIA1',
              'X_LIA1_corr', 'LIA1_DMM', 'LIA1corr']
    Gz = abs(hedict['G_zero'])
    if 1E-4 <= Gz < 1E-2:
        fac = 1E3
        Gaxis = '$dI/dV \ [mG/G_{0}]$'
    elif 1E-7 <= Gz < 1E-4:
        fac = 1E6
        Gaxis = '$dI/dV \ [ \u00b5G/G_{0}]$'
    elif Gz < 1E-7:
        fac = 1E9
        Gaxis = '$dI/dV \ [nG/G_{0}]$'
    else:
        fac = 1
        Gaxis = '$dI/dV \ [G/G_{0}]$'
    for m, n in enumerate(Dictionaries):
        for i, j in enumerate(G_keys):
            if j in n.keys():
                n[j] = n[j]*fac
    return Dictionaries, Gaxis, fac


def save_dictionary(dictio, path, name):
    #mpal = makepath_local(paths, name, excu, usename = 'yes')
    for i, j in enumerate(dictio):
        fn = path + '/' + name[i] + '.json'
        with open(fn, "w") as fi:
            json.dump(j, fi)
        fi.close()


def copy_files(src_lst, dst_):
    for i, j in enumerate(src_lst):
        name = os.path.split(j)
        sht.copy2(j, dst_ + '/' + name[1])


def tolist(dcts):
    n_lst = []
    for i, j in enumerate(dcts):
        n_dc = {}
        for l, m in enumerate(j.keys()):
            n_dc[m] = list(j[m])
        n_lst.append(n_dc)
    return n_lst


def smooth_param(x, sm_degree):
    """ degree of smoothin can be 'rough', 'middle', or 'fine'.
    """
    lx = len(x)
    if lx < 150:
        lx = lx + (150 - lx)
    if sm_degree == 'rough':
        win = round(lx/120)
        if win % 2 == 0:
            win = win+1
        polyorder = 2
        deriv = 0
        delta = 0
        last = 'interp'
    elif sm_degree == 'middle':
        win = round(lx/40)
        if win % 2 == 0:
            win = win+1
        polyorder = 4
        deriv = 0
        delta = 0
        last = 'interp'
    elif sm_degree == 'fine':
        win = round(lx/15)
        if win % 2 == 0:
            win = win+1
        polyorder = 6
        deriv = 0
        delta = 0
        last = 'interp'
    return [win, polyorder, deriv, delta, last]


def smoothing(y, filterparams):  # =[21,0,0,0.5,'interp']):
    """y can be a List of arrays"""
    window_sm = filterparams[0]
    polyorder = filterparams[1]
    deriv = filterparams[2]
    delta_sm = filterparams[3]
    mode = filterparams[4]
    Smi = signal.savgol_filter(y, window_sm,
                               polyorder, deriv, delta_sm, mode=mode)
    return Smi


def numeric(dct, key, gfac, smd='rough'):
    'smd -- > smooth degree'
    V = dct[key[0]]
    I = smoothing(dct[key[1]], smooth_param(V, smd))
    LIA1 = smoothing(dct[key[2]], smooth_param(V, smd))
    LIA2 = smoothing(dct[key[3]], smooth_param(V, smd))
    G0 = spc.physical_constants['conductance quantum'][0]
    #Glin = np.mean(np.isclose(I, 0, atol=1))
    dIdV = np.gradient(I/1E6, V/1E3)/G0*gfac
    dIdV_ = integ.cumtrapz(LIA2, V/1E3)*gfac
    dGdV = np.gradient(LIA1/gfac, V/1E3)
    dIdVsm = smoothing(dIdV, smooth_param(dIdV, smd))
    dGdVsm = smoothing(dGdV, smooth_param(dGdV, smd))
    I = integ.cumtrapz(LIA1*G0, V/1E3)*1E6/gfac
    return [I, dIdVsm, dIdV_, dGdVsm]


def keymake(dct, corr='yes'):
    if corr == 'yes':
        if 'LIA_1' or 'LIA1' or 'LIA_1_DMM' in dct.keys():
            keys = ['Bias_DMM', 'I(DMM)', 'LIA1corr', 'LIA2_DMM', 'XY']
        else:
            keys = ['Bias_DMM', 'I(DMM)', 'X_LIA1corr', 'LIA2_DMM', 'XY']
    else:
        if 'LIA_1' or 'LIA1' or 'LIA_1_DMM' in dct.keys():
            keys = ['Bias_DMM', 'I(DMM)', 'LIA1_DMM', 'LIA2_DMM', 'XY']
        else:
            keys = ['Bias_DMM', 'I(DMM)', 'X_LIA1', 'LIA2_DMM', 'XY']
    return keys


def sweepdata(idx, Headers, file):
    he = Headers[0]
    sweepdat = [idx, he['conductance'], he['Motor_Pos'], he['Data_Points'],
                he['O/C_Cyles'], file]
    return sweepdat

# figures =====================================================================


def figcheck(dct, key, vboc, Goc, ylab, leg):
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.scatter(vboc, Goc, s=5, c='r')
    ax.plot(dct[key[0]], dct[key[1]])
    ax.set_ylabel(ylab, fontsize=14)
    ax.set_xlabel('Bias [mV]', fontsize=14)
    ax.legend([leg])
    f.tight_layout()
    return f


def fig_stack(dct, leg, keys, xlab, ylab):
    """ signal 'iets' or 'dG/dV' """
    plt.close('all')
    #I,LIA1,LIA2, Iten = y[0], y[1], y[2], y[3]
    f1 = plt.figure(1, figsize=(6, 11))
    ax4 = f1.add_subplot(414)
    ax1 = f1.add_subplot(411, sharex=ax4)
    ax2 = f1.add_subplot(412, sharex=ax4)
    ax3 = f1.add_subplot(413, sharex=ax4)
    #ax4 = ax3.twinx()
    f1.tight_layout()
    ax1.plot(dct[keys[0]], dct[keys[1]], c='k')
    ax1.set_ylabel('current [\u00b5A]', fontsize=14)
    ax1.tick_params(axis='both', labelsize=12.0)
    ax2.plot(dct[keys[0]], dct[keys[2]], c='k')
    ax2.set_ylabel(ylab, fontsize=14)
    ax2.tick_params(axis='both', labelsize=12.0)
    ax3.plot(dct[keys[0]], dct[keys[3]], c='k')  # , label='dG/dV [G0]')
    ax3.set_ylabel('$dG/dV \ [(G/G_{0})*V^{-1}]$', fontsize=14)
    ax4.plot(dct[keys[0]], dct[keys[4]], c='k', alpha=1)
    ax4.set_ylabel('$Intensity \ [V^{-1}]$', fontsize=14)
    ax3.set_xlabel('Bias [mV]', fontsize=14)
    ax3.tick_params(axis='both', labelsize=12.0)
    ax3.legend([leg])
    f1.tight_layout()
    return f1


def fig_ad_sub(V, y, fig, ax, plotpos, x_label, y_label, cl=[]):
    """ads a subplot to figure 'fif' at postition 'plotpos'
    plotpos is given as triple 'eg. 111' or [111] with the axis name ax
    cl:label of the curve if not then nothin happens
    """
    ax = fig.add_subplot(plotpos)
    c = ['k', 'b', 'm', 'c']
    a = [0.8, 1, 1, 1]
    for i in range(len(V)):
        if len(cl) > 0:
            ax.plot(V[i], y[i], color=c[i], alpha=a[i], label=cl[i])
            handles, labels = ax.get_legend_handles_labels()
            plt.legend()
        else:
            ax.plot(V[i], y[i], color=c[i], alpha=a[i])
    ax.set_xlabel(x_label, fontsize=16)
    ax.set_ylabel(y_label, fontsize=16)
    ax.tick_params(axis='both', labelsize=12.0)
    return ax


def fig_up_down(dicts, gylab, key):
    """ plots is array with the number of plots and what should be plotted.
    it can be ['I','dI/dV', 'dG/dV'] (iets or dG/dv) or just one or two of
    them. The last entry in the list is the lower plot. cl is the label
    of the curves if there are more than one in a plot.
    """
    #y,yn,yp = y_all
    #I,LIA1,LIA2, iets = y[0], y[1], y[2], y[3]
    #In,LIA1n,LIA2n,ietsn = yn[0], yn[1], yn[2], yn[3]
    # Ip,LIA1p,LIA2p,ietsp = yp[0], yp[1], yp[2], yp[3]list(thet)
    #V, Vn, Vp = V_all

    Is = [i[key[1]] for i in dicts]
    LIA1s = [i[key[2]] for i in dicts]
    LIA2s = [i[key[3]] for i in dicts]
    ietss = [i[key[4]] for i in dicts]
    Vs = [i[key[0]] for i in dicts]

    f_ud = plt.figure(1, figsize=(5.6, 9))
    fig_ad_sub(Vs, Is, f_ud, 'ax1', 411, 'Bias [mV]', 'current [\u00b5A]')
    fig_ad_sub(Vs, LIA1s, f_ud, 'ax2', 412, 'Bias [mV]', gylab)
    fig_ad_sub(Vs, ietss, f_ud, 'ax4', 414, 'Bias [mV]', 'IETS [1/V]')
    fig_ad_sub(Vs, LIA2s, f_ud, 'ax3', 413, 'Bias [mV]',
               '$dG/dV \ [(G/G_{0})*V^{-1}]$')
    f_ud.tight_layout()
    return f_ud


def fig_deriv(dct, keys, num, ylab, gz, legend='yes', lim=7):
    """if an integration is calculated"""
    V, I, LIA1, LIA2 = dct[keys[0]], dct[keys[1]], dct[keys[2]], dct[keys[3]]
    I_in, dIdVd, dIdVI_, dGdV = num
    le = round(len(V)/2)
    I_int = I_in - I_in[le]
    dIdVI = dIdVI_ + gz
    plt.close('all')
    f2 = plt.figure(2, figsize=(6, 8.5))
    ax3 = f2.add_subplot(313)
    ax1 = f2.add_subplot(311, sharex=ax3)
    ax2 = f2.add_subplot(312, sharex=ax1)
# =============================================================================
#     ax11 = ax1.twinx()
#     ax22 = ax2.twinx()
#     ax31 = ax3.twinx()
# =============================================================================
    ax1.plot(V, I, c='k', label='measurement')
    ax1.plot(V[:-1], I_int, '--', c='c',  label='integration')
    ax1.set_ylabel('current [\u00b5A]', fontsize=16)
    ax1.tick_params(axis='both', labelsize=12.0)

    ax2.plot(V, LIA1, c='k', label='measurement')
    ax2.plot(V[:-1], dIdVI, '--', c='c', label='integration')
    ax2.plot(V[lim:-(lim)], dIdVd[lim:-lim], '--', c='m', label='derivative')
    ax2.set_ylabel(ylab, fontsize=16)
    ax2.tick_params(axis='both', labelsize=12.0)

    ax3.plot(V, LIA2, c='k')
    ax3.plot(V[lim:-(lim)], dGdV[lim:-lim], '--', c='c', label='derivative')
    ax3.set_ylabel('$dG/dV \ [(G/G_{0})*V^{-1}]$', fontsize=16)
    ax3.set_xlabel('Bias [mV]', fontsize=16)
    ax3.tick_params(axis='both', labelsize=12.0)

    ax1.legend()
    ax3.legend()
    ax2.legend()

    f2.tight_layout()

    return f2


def fig_sing(dct, key, labels, leg=None):

    f = plt.figure()
    ax = f.add_subplot(111)
    ax.plot(dct[key[0]], dct[key[1]], color='grey', alpha=0.85)
    ax.scatter(dct[key[0]], dct[key[1]], s=1, c='navy', alpha=1)
    ax.set_ylabel(labels[1], fontsize=16)
    ax.set_xlabel(labels[0], fontsize=16)
    ax.legend([leg])
    f.tight_layout()
    return f


def fig_sing_updown(dct, key, labels):
    plt.close('all')
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.plot(dct[0][key[0]], dct[0][key[1]], label='full sweep')
    ax.plot(dct[1][key[0]], dct[1][key[1]], label='first sweep')
    ax.plot(dct[2][key[0]], dct[2][key[1]], label='second sweep')
    ax.set_ylabel(labels[1], fontsize=16)
    ax.set_xlabel(labels[0], fontsize=16)
    ax.legend()
    f.tight_layout()
    return f


def fig_sing_lim(dct, key, labels, lim, leg=None):
    V = dct[key[0]]
    #V1 = dct[1][key[0]]
    #V2 = dct[2][key[0]]
    mask = (V >= lim[0]) & (V <= lim[1])
    #mask1 = (V1 >= lim[0]) & (V1 <= lim[1])
    #mask2 = (V2 >= lim[0]) & (V2 <= lim[1])
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.scatter(dct[key[0]][mask], dct[key[1]][mask], s=1, c='r')
    ax.plot(dct[key[0]][mask], dct[key[1]][mask], color='k')
    ax.set_ylabel(labels[1], fontsize=16)
    ax.set_xlabel(labels[0], fontsize=16)
    ax.xaxis.set_major_locator(MultipleLocator(50))
    ax.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    plt.grid(True, which='major', axis='x', color='r', linewidth=0.5)
    plt.grid(True, which='minor', axis='x')
    ax.legend([leg])
    f.tight_layout()
    return f


def fig_scat_updown(dct, key, labels):
    plt.close('all')
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.scatter(dct[0][key[0]], dct[0][key[1]], label='full sweep', s=7)
    ax.scatter(dct[1][key[0]], dct[1][key[1]], label='first sweep', s=4)
    ax.scatter(dct[2][key[0]], dct[2][key[1]], label='second sweep', s=4)
    ax.set_ylabel(labels[1], fontsize=16)
    ax.set_xlabel(labels[0], fontsize=16)
    ax.legend()
    f.tight_layout()
    return f


def G_colorjet(x, y, param, rang, title, design='bone'):
    """ signal 'iets' or 'dGdV', I or didv """
    plt.close('all')
    norm = mpl.colors.Normalize(vmin=4, vmax=max(param)*1.4)
    f1 = plt.figure(3, figsize=[9.6, 4.8])
    ax1 = f1.add_subplot(111)
    plt.ylim(bottom=min(y)-min(y)/2)
    plt.ylim(top=max(y) + max(y)/2)
    ax1.set_yscale('log')
    f1.tight_layout()
    plot = ax1.scatter(x, y, c=param, cmap=design,
                       norm=norm, s=0.5, alpha=0.95)
    # ax1.pcolor([V,param])
    cbar = f1.colorbar(plot)
    cbar.set_label('Temperature [K]', size=16)
    ax1.set_title(title, size=16)
    #f1.subtitle('OV Bias Conductivity < 0.05 G0',fontsize=16)
    ax1.set_ylabel('$G \/\ [G/G_0]$', fontsize=18)
    ax1.set_xlabel('Motorposition [a.u.]', fontsize=18)
    ax1.tick_params(axis='both', labelsize=16.0)
    if MposOC[0] > MposOC[len(MposOC)-1]:
        ax.invert_xaxis()
    #ax1.legend(['N = ' + str(len(y))+ ' points'])
    f1.tight_layout()
    return f1


def G_colorjet_IV(x, y, swpos, param, rang, title, design='bone'):
    """ signal 'iets' or 'dGdV', I or didv """
    #Mpos_OC = dct['MposOC']
    plt.close('all')
    norm = mpl.colors.Normalize(vmin=4, vmax=max(param)*1.4)
    f1 = plt.figure(3, figsize=[9.6, 4.8])
    ax1 = f1.add_subplot(111)
    plt.ylim(bottom=min(y)-min(y)/2)
    plt.ylim(top=max(y) + max(y)/2)
    ax1.set_yscale('log')
    f1.tight_layout()

    plot = ax1.scatter(x, y, c=param, cmap=design,
                       norm=norm, s=0.5, alpha=0.95)

    # ax1.pcolor([V,param])
    cbar = f1.colorbar(plot)
    cbar.set_label('Temperature [K]', size=18)
    ax1.set_title(title, size=16)
    #f1.subtitle('OV Bias Conductivity < 0.05 G0',fontsize=16)
    ax1.scatter((swpos[0]/1E+7), swpos[1], color='r', s=1.5)
    ax1.set_ylabel('$G \/\ [G/G_0]$', fontsize=20)
    ax1.set_xlabel('Motorposition [a.u.]', fontsize=20)
    ax1.tick_params(axis='both', labelsize=18.0)
    if MposOC[0] > MposOC[len(MposOC)-1]:
        ax1.invert_xaxis()
    #ax1.legend(['N = ' + str(len(y))+ ' points'])
    f1.tight_layout()
    return f1


def fig_compare(dct, allsweepdat):
    plt.close('all')
    fig, ax = plt.subplots(figsize=(9.6, 4.8))
    Gp, = ax.plot(dct['MposOC'], dct['Gval'],  color='b')
    Gpd, = ax.plot(dct['MposOC'], dct['LIAG'], color='m')
    if min(dct['Gval']) > min(dct['LIAG']):
        plt.ylim(bottom=min(dct['LIAG'])-min(dct['LIAG'])/2,)
    else:
        plt.ylim(bottom=min(dct['Gval'])-min(dct['Gval'])/2, )
    if max(dct['Gval']) > max(dct['LIAG']):
        plt.ylim(top=max(dct['Gval'])+max(dct['Gval'])/2)
    else:
        plt.ylim(top=max(dct['LIAG'])+max(dct['LIAG'])/2)
    ax.set_yscale('log')
    plt.grid(True, which='major', axis='y', color='g', linewidth=0.2)
    ax.set_xlabel('Motorposition', fontsize=18)
    ax.set_ylabel('conductance [G0]', fontsize=18)
    ax.tick_params(axis='both', labelsize=16.0)
    ax2 = ax.twinx()
    ax2.scatter(dct['MposOC'], dct['vel'], s=0.1, color='k', alpha=0.15)
    plt.legend([Gp, Gpd], ['static G', 'differential G'], fontsize=14,
               loc='best')
    ax2.tick_params(axis='y', labelsize=10.0)

    # if abs(MposOC[0]) > abs(MposOC[-1]):
    if dct['MposOC'][0] > dct['MposOC'][-1]:
        ax.invert_xaxis()
    swpos = []
    for i in range(len(allsweepdat)):
        if i > 0:
            if i == allsweepdat[i][4]:
                swpos.append([allsweepdat[i][2], allsweepdat[i][1]])
    swpos = np.transpose(swpos)
    ax.scatter((swpos[0]/1E+7), swpos[1], color='r')
    plt.tight_layout()
    return fig