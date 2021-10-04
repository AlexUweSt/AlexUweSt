# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 15:56:21 2021

@author: strobe46
"""

import datetime as dt
import os
import scipy.constants as spc
import scipy.integrate as integ
import matplotlib.pyplot as plt
import numpy as np
import json
from scipy.fft import fft, ifft
from scipy.signal import blackman, blackmanharris, hann, tukey
import peakdetect as pedet
import math
from scipy import signal
from lmfit import Model
from lmfit.models import LorentzianModel, QuadraticModel, GaussianModel,\
    LinearModel, ConstantModel, SplitLorentzianModel, VoigtModel,\
        BreitWignerModel, SkewedGaussianModel,DampedHarmonicOscillatorModel,\
            SkewedVoigtModel,ThermalDistributionModel, DoniachModel,\
                LognormalModel#, FanoModel, ExponentialGaussianModel
from scipy import interpolate
##########################################################################

def load_full_analyze(fpath, name, cutar,xfac, Acamp, ti):
    fpath = fpath + '/processed_arrays_files'
    G0 = spc.physical_constants['conductance quantum'][0]
    MEAS = np.loadtxt(fpath + '/fullsweep.txt' ,dtype=float, delimiter = ' ' )
    MEAS_T=np.transpose(MEAS)
    Header=np.genfromtxt(fpath + '/Header.txt' ,dtype='str',max_rows=6)
    Header = np.char.replace(Header, ',', '.')
    if Acamp == 'auto':
        #ACamp = int(Header[3][6])/1000   #in volts
        ACamp = float(Header[3][6])   #in volts
    else:
            ACamp = Acamp
    if ti == 'yes':
        tistamp = np.load(fpath + '/time.npy', allow_pickle=True)
        #times = [tistamp[i].time() for i in tistamp]
        ti = [i.timestamp() for i in tistamp]
        sec = [i-min(ti) for i in ti]
        sec = np.array(sec)
    else:
        sec = [0,0]
    # =============================================================================
    # Pkdiff = (np.sqrt((5.4*spc.Boltzmann*T_insitu/spc.e)**2 +\
    #                  (1.7*ACamp)**2))*1.00 #sqrt(x) is a random factor..
    # =============================================================================
    param = Header[1].astype(float)

    R_pre  = float(Header[1][0]) #in Ohms
    R_pre_G = (1/R_pre)/G0
    #Vad = np.loadtxt(fpath+'/V_IETS.txt')
    Vad = MEAS_T[8][cutar:]
    Vspt = MEAS_T[0][cutar:] #V_setpoint
    vdrop = MEAS_T[5][cutar:]
    us = (Vad-(vdrop*param[7]*R_pre/2)) # voltagedrop Usample
    I=Vad/(us/(vdrop*param[7])) * (1E6)
    Imin = min(abs(I))
    Imin_ind, = np.where(abs(I)==Imin)[0]
    if float(Header[1,8]) < 0:
        fac1 = 1
    else:
        fac1 = -1

    if float(Header[3,3]) < 0:
        fac2 = 1
    else:
        fac2 = -1

    LIA1_gain = param[3]/(param[4]*10) * fac1
    LIA1_a = MEAS_T[6][cutar:]*param[7]/(ACamp*G0)*LIA1_gain
    LIA1 = 1/((1/LIA1_a) - (1/R_pre_G))
    LIA2_gain = param[5]/(param[6]*10) * fac2
    LIA2 = MEAS_T[7][cutar:]*param[7]/(ACamp*ACamp*G0)*LIA2_gain
    Iten = LIA2/LIA1
    y = [I,LIA1,LIA2,Iten]

    #load up and downsweep:
    if 'neghalfsweep.txt' in os.listdir(fpath):
        neg = np.loadtxt(fpath + '/neghalfsweep.txt')#,dtype=float, delimiter = ' \t' )
    elif '1st_sweep.txt' in os.listdir(fpath):
        neg =  np.loadtxt(fpath + '/1st_sweep.txt')
    #neg=np.transpose(neg)
    Vn = neg[8]
    Vns = (Vn-(neg[5]*param[7]*R_pre)) # voltagedrop
    In = Vn/((Vn-(neg[5]*param[7]*R_pre))/(neg[5]*param[7])) * (1E6)
    LIA1n_a = neg[6]*param[7]/(ACamp*G0)*LIA1_gain
    LIA1n = LIA1 = 1/((1/LIA1n_a) - (1/R_pre_G))
    LIA2n = neg[7]*param[7]/(ACamp*ACamp*G0)*LIA2_gain
    Itenn = LIA2n/LIA1n
    yn = [In,LIA1n,LIA2n,Itenn]
    #oc_bias = round(Vn[0],2)*1000
    if 'poshalfsweep.txt' in os.listdir(fpath):
        pos = np.loadtxt(fpath + '/poshalfsweep.txt')# ,dtype=float, delimiter = ' \t' )
    elif '2nd_sweep.txt' in os.listdir(fpath):
        pos = np.loadtxt(fpath + '/2nd_sweep.txt')
    #pos=np.transpose(pos)
    Vp = pos[8]
    Vps = (Vp-(pos[5]*param[7]*R_pre))
    Ip = Vp/((Vp-(pos[5]*param[7]*R_pre))/(pos[5]*param[7])) * (1E6)
    LIA1p_a = pos[6]*param[7]/(ACamp*G0)*LIA1_gain
    LIA1p = 1/((1/LIA1p_a) - (1/R_pre_G))
    LIA2p = pos[7]*param[7]/(ACamp*ACamp*G0)*LIA2_gain
    Itenp = LIA2p/LIA1p
    yp = [Ip,LIA1p,LIA2p,Itenp]
    Vall = [us*xfac,Vns*xfac,Vps*xfac]
    yall = [y,yn,yp]
    
    return Vall, yall, Header, sec

def load_processed(fp, name, cutar, ti='no', onefile = 'true'):
    fpat = fp + '/processed_arrays_files'
    if onefile =='false':
        with open(fpat + '/fullsweep.json', 'r') as file:
            full=file.read()
        Meas = json.loads(full)#,allow_pickle=True)
        
        with open(fpat + '/fullsweep.json', 'r') as file:
            fi=file.read()
        first = json.loads(fi)#,allow_pickle=True)
        
        with open(fpat + '/fullsweep.json', 'r') as file:
            se=file.read()
        secon = json.loads(se)#,allow_pickle=True)
        
        with open(fpat + '/header_dict.json', 'r') as file:
            he=file.read()
        Header = json.loads(he)#,allow_pickle=True)    
        
    elif onefile == 'true':
        with open(fpat + '/meas_files.json', 'r') as file:
            full=file.read()
        afi = json.loads(full)#,allow_pickle=True)     # afi -> allfiles
    Meas, first, secon, Header = afi['full'], afi['first'], afi['second'],\
        afi['header']
    
    if ti == 'yes':
        tistamp = Meas['time']
        #times = [tistamp[i].time() for i in tistamp]
        ti = [i.timestamp() for i in tistamp]
        sec = [i-min(ti) for i in ti]
        sec = np.array(sec)
    else:
        sec = [0,0]
    
    #create Arrays:
    if 'LIA1corr' in Meas.keys():
        lia1 = 'LIA1corr'
    else: 
        lia1 = 'X_LIA1corr'
        
    y = [np.array(Meas['I(DMM)'][cutar:]), np.array(Meas[lia1][cutar:]),\
         np.array(Meas['LIA2_DMM'][cutar:]), np.array(Meas['XY'][cutar:])]
    y1 = [np.array(first['I(DMM)'][cutar:]), np.array(first[lia1][cutar:]),\
          np.array(first['LIA2_DMM'][cutar:]), np.array(first['XY'][cutar:])] 
    y2 = [np.array(secon['I(DMM)'][cutar:]), np.array(secon[lia1][cutar:]),\
          np.array(secon['LIA2_DMM'][cutar:]), np.array(secon['XY'][cutar:])]
    Vall = [np.array(Meas['Bias_DMM'][cutar:]), np.array(first['Bias_DMM'][cutar:]),\
            np.array(secon['Bias_DMM'][cutar:])]
    yall = [y, y1, y2]
    
    return Vall, yall, Header, sec

def changename(oldname, ad):
    newad = ad.replace('/', '_')
    newname = newad+ oldname
    return newname

def Parameters(Header, dictio):
    if type(Header) != dict:
        Legend = Header[1,1] +' G0' + ' @ '+ Header[3][6]+ '°C'
        T_insitu = float(Header[3][4]) # int(Header[3][4])
        G_OC = Header[1][1]
        Lockamp = Header[3][6]
        ad_header_to_dict(Header,'Header',dictio)
    else: 
        T_insitu = float(Header['Temperature']) # int(Header[3][4])
        G_OC = Header['G_zero']
        Legend = str(G_OC) +' G0' + ' @ '+ str(T_insitu)+ '°C'
        Lockamp = Header['LockAmpl']
        dictio['Header'] = Header
    diff = (np.sqrt((5.4*spc.Boltzmann*T_insitu/spc.e)**2 +\
             (1.7*float(Lockamp)/1000)**2))*np.sqrt(1.00) #sqrt(x) is a random factor..
    return Legend, T_insitu, G_OC, diff

def make_plots_path(fp, name = 'analyzer_plots'):
    if os.path.exists(fp+'/' + name) == False:
        os.mkdir(fp + '/' + name)
    return fp + '/' + name

def makepath_local(paths, name, excup, usename = 'no'):
    """ excup is atrributes, extrema and curves path, where entry [0] is
    the attributes, entry [1] is extrema and entry [2] is the curves path
    """
    if type(excup) == list:
        exp = excup[0]
        cup = excup[1]
        atp = excup[2]
    elif type(excup) == str:
        exp = cup = atp = excup
    cupa, expa, atpa = [], [], []
    cupp, expp, atpp = [], [], []
    for i, j in enumerate(paths):
        if usename == 'yes':
            name0 = j + '/' + exp + '/' + name
            name1 = j + '/' + cup + '/' + name
            name2 = j + '/' + atp + '/' + name
            expa.append(name0)
            cupa.append(name1)
            atpa.append(name2)
        elif usename == 'no':
            name0 = j+ '/' + exp
            name1 = j+ '/' + cup
            name2 = j+ '/' + atp
            expa.append(name0)
            cupa.append(name1)
            atpa.append(name2)
        name00 = j+ '/' + exp
        name10 = j+ '/' + cup
        name20 = j+ '/' + atp
        expp.append(name00)
        cupp.append(name10)
        atpp.append(name20)
    for m,n in enumerate(atpp):
        if os.path.exists(n) == False:
            #pa = Path(n)
            os.makedirs(n)
        elif os.path.exists(n) == True:
                pass
    for m,n in enumerate(expp):
        if os.path.exists(n) == False:
            #pa = Path(n)
            os.makedirs(n)
        elif os.path.exists(n) == True:
                pass
    for m,n in enumerate(cupp):
        if os.path.exists(n) == False:
            #pa = Path(n)
            os.makedirs(n)
        elif os.path.exists(n) == True:
                pass
    return expa, cupa, atpa

def save_dictionaries(dicts, paths, name, excup):
    mpal = makepath_local(paths, name, excup, usename = 'yes')
    for i, j in enumerate(mpal[0]):
        fn = j + '.json'
        with open(fn, "w") as fi:
            json.dump(dicts[0], fi)
    fi.close()
    for i, j in enumerate(mpal[1]):
        fn = j + '.json'
        with open(fn, "w") as fi:
            json.dump(dicts[1], fi)
    fi.close()
    for i, j in enumerate(mpal[2]):
        fn = j + '.json'
        with open(fn, "w") as fi:
            json.dump(dicts[2], fi)
    fi.close()


def ad_header_to_dict(header,name,dictio):
    di = {}
    for i,j in enumerate(header[0]):
        di[j] = float(header[1][i])
        di[header[2][i]] = float(header[3][i])
    dictio[name] = di

def titofr(arr):  # time to frequency ..
    freq = []
    dta = []
    for i in range(len(arr[:-1])):
        delta = arr[i+1]-arr[i]
        fr = 1/delta
        dta.append(delta)
        freq.append(fr)
    dta = np.array(dta)
    freq = np.array(freq)
    return freq, dta

def ycutter(y, cut):
    if type(y)==list:
        yc = []
        for i in range(len(y)):
            yc.append(y[i][cut])
    elif type(y)==np.ndarray:
        yc = y[cut]
    return yc

def G_y_label(G_fac):
    if G_fac == 1E3:
        gylab = '$dI/dV \ [mG/G_{0}]$'# '$dI/dV \ [G/G_{0}]$' # #
    elif G_fac == 1E0:
        gylab = '$dI/dV \ [G/G_{0}]$'
    elif G_fac == 1E6:
        gylab = '$dI/dV \ [ \u00b5G/G_{0}]$'
    elif G_fac == 1E9:
        gylab = '$dI/dV \ [nG/G_{0}]$'
    return gylab
    

def subs_lin_bkg(x,y,c,fac):
    dx = abs(max(x)-min(x))
    dy = abs(max(y) - min(y))
    a = (dy/dx)*fac
    ar = np.linspace(min(x), max(x), len(x))
    #fc = a*r
    return a*ar+c, [a,c]

def IV_cutter(t,y_2darr, cutout):
    y = y_2darr
    if type(t) == np.ndarray:
        mI = (t >= cutout[0]) & (t <= cutout[1]) #mask
        yc = ycutter(smooth(y,smooth_param(t,'rough')), mI)#[y_arr]
        tc = t[mI]
    elif type(t) == list:
        mI = (t[0] >= cutout[0]) & (t[0] <= cutout[1]) # --> mask1
        yc = ycutter(smooth(y,smooth_param(t[0],'rough')), mI)#[y_arr]
        tc = []
        for i in range(len(t)):
            tc.append(t[i][mI])
    return tc,yc

def fourier(t, y,cut, wfct='hann', f_win=[5E-3,0.035],ampf=1E3):
    """ampf is the amplitude amplification factor
    f_win is the considered frequency window, maximum should depend on
    the Nyquist–Shannon sampling theorem
    wfct is the window funtion applied on the y axis array to exlude effects
    on the edges --> lackman, blackmanharris, hann, tukey are the most
    reasonable ones for this application (I think..)"""
    tc, yc = IV_cutter(t,y, cut, 1)
    if type(yc) == list:
        yc = yc[1]
    ###### window functions
    #bmw = tukey(zN)
    bmw = hann(int(len(yc)))
    #bmw = blackman(zN)
    # fourier transform
    ffty = fft(yc)
    #r_fft = np.abs(ffty)
    prod = yc*bmw
    w_fft = fft(prod)

    #prepare the frequeny vector
    freq = titofr(t)[0]
    #avdT = round(np.average(ddt))
    #avfr = np.average(freq)
    nq = min(freq)/2
    N = int(len(yc)/2+1)
    fvec = np.linspace(0,nq,N)
    mII = (fvec >= f_win[0]) & (fvec <= f_win[1])
    fvec = fvec[mII]  #*1E3
    refftw = ampf*(2/int(len(yc)))*np.abs( w_fft[0:N])[mII]
    refft = ampf*(2/int(len(yc)))*np.abs( ffty[0:N])[mII]
    maxfr = fvec[np.where(refftw == max(refftw))]
    return refft, refftw, fvec, tc, yc, prod, maxfr

def integrator(x,y,facx=1000, facy = 1, c='dgdv'):
    """
    ..

    """
    #G0 = spc.physical_constants['conductance quantum'][0]
    xa = np.argsort(x)
    x_ = np.array(x)[xa]/facx
    y_ = np.array(y)[xa]/facy
    dct = {}
    [Vp, yp], [Vn_, yn_] = Reflection(x_,y_)[1:]
    Vn, yn = [abs(i) for i in Vn_][::-1], yn_[::-1]
    # = Vp_[::-1], yp_[::-1]
    Vra = abs(min(x_)) + abs(max(x_))#Bias range all
    Vrpn = abs(max(Vn)) + abs(max(Vp))
    Vrp = abs(max(Vp))
    Vrn = abs(max(Vn))
    int_a = integ.simps(y_,x_,even='avg')#/G0
    int_p = integ.simps(yp,Vp,even='avg')#/G0
    int_n = integ.simps(yn,Vn,even='avg')#/G0
    if c == 'dgdv':
        dct['dgdv_neg+pos'] = int_p + int_n, (int_p + int_n)/Vrpn
        dct['dgdv_whole_range'] = int_a, int_a/Vra
        dct['dgdv_pos_arm'] = int_p, int_p/Vrp
        dct['dgdv_neg_arm'] = int_n, int_n/Vrn
        dct['dgdv_ratio_pn']  = (int_p/Vrp)/(int_n/Vrn)
        dct['dgdv_factors_x_Y'] = facx,facy
    elif c == 'iets':
        dct['iets_neg+pos'] = int_p + int_n, (int_p + int_n)/Vrpn
        dct['iets_whole_range'] = int_a, int_a/Vra
        dct['iets_pos_arm'] = int_p, int_p/Vrp
        dct['iets_neg_arm'] = int_n, int_n/Vrn
        dct['iets_ratio_pn']  = (int_p/Vrp)/(int_n/Vrn)
        dct['iets_factors_x_Y'] = facx,facy
    elif c == 'didv':
        dct['didv_neg+pos'] = int_p + int_n, (int_p + int_n)/Vrpn
        dct['didv_whole_range_'] = int_a, int_a/(Vra)#*G0)
        dct['didv_pos_arm_'] = int_p, int_p/(Vrp)#*G0)
        dct['didv_neg_arm_'] = int_n, int_n/(Vrn)#*G0)
        dct['didv_ratio_pn_']  = (int_p/Vrp)/(int_n/Vrn)
        dct['didv_factors_x_Y'] = facx,facy
    dct['ranges_n_p_all'] = Vrn, Vrp, Vra
    return dct

def symm_attr(x,y,curve='dIdV', facx=1E3,facy=1.):
    y = y/facy
    x = x/facx
    dct = {}
    px, py = Reflection(x,y)[1]
    nx, ny = Reflection(x,y)[2]
    #dct = np.where(ny == max(ny))[0][0]
    if curve == 'dIdV':
        dct[curve + 'multply_fac_x_y'] = facx, facy
        dct[curve+'max_global']= mg =(x[np.where(y == max(y))[0][0]],max(y))
        dct[curve+'max_neg'] = mn =(nx[np.where(ny == max(ny))[0][0]],max(ny))
        dct[curve+'max_pos'] = mp =(px[np.where(py == max(py))[0][0]],max(py))
        dct[curve+'min_global'] = mig =(x[np.where(y == min(y))[0][0]],min(y))
        dct[curve+'min_neg'] =minn=(nx[np.where(ny == min(ny))[0][0]],min(ny))
        dct[curve+'min_pos'] = mip=(px[np.where(py == min(py))[0][0]],min(py))
        dct[curve+'onoff_global'] = mg[1]/mig[1]
        dct[curve+'onoff_pos'] = mp[1]/mip[1]
        dct[curve+'onoff_neg'] = mn[1]/minn[1]
    elif curve == 'I':
        dct[curve + 'multply_fac_x_y'] = facx, facy
        dct[curve+'max_neg'] = mn =(nx[np.where(ny == max(ny))[0][0]],max(ny))
        dct[curve+'max_pos'] = mp =(px[np.where(py == max(py))[0][0]],max(py))
        dct[curve+'min_neg'] =minn=(nx[np.where(ny == min(ny))[0][0]],min(ny))
        dct[curve+'min_pos'] = mip=(px[np.where(py == min(py))[0][0]],min(py))
    elif curve == 'dGdV' or 'IETS':
        dct[curve + 'multply_fac_x_y'] = facx, facy
        dct[curve+'max_global']= mg =(x[np.where(y == max(y))[0][0]],max(y))
        dct[curve+'max_neg'] = mn =(nx[np.where(ny == max(ny))[0][0]],max(ny))
        dct[curve+'max_pos'] = mp =(px[np.where(py == max(py))[0][0]],max(py))
        dct[curve+'min_global'] = mig =(x[np.where(y == min(y))[0][0]],min(y))
        dct[curve+'min_neg'] =minn=(nx[np.where(ny == min(ny))[0][0]],min(ny))
        dct[curve+'min_pos'] = mip=(px[np.where(py == min(py))[0][0]],min(py))
    return dct

def prep_Coul_(pd, G_di, pdi):
    """
    Parameters
    ----------
    pd : peaks, dips (list or numpy array)
    G_di : dictionary
        DESCRIPTION.
    pdi : list of 2 lists
        indeces to pick out the interesting peaks and dips.

    Returns
    -------
    di : TYPE
        DESCRIPTION.

    """
    peaks_, dips_ = pd
    xpeaks = peaks_[0][pdi[0]] #[peaks_[0][0],peaks_[0][-1]]
    ypeaks = peaks_[1][pdi[0]] #[peaks_[1][0],peaks_[1][-1]]
    xdips = dips_[0][pdi[1]] #[peaks_[0][0],peaks_[0][-1]]
    ydips = dips_[1][pdi[1]] #[peaks_[1][0],peaks_[1][-1]]
    ylev = np.average(ypeaks)
    peaks = np.vstack((xpeaks,ypeaks)) # should be 2 points
    dips = np.vstack((xdips,ydips)) # should be a single point
    exx = np.concatenate((peaks[0],dips[0])) # extrema x axis
    exy = np.concatenate((peaks[1],dips[1])) # extrema y-axis
    ex = np.vstack((exx,exy)).tolist() # extrema
    p0, p1 = round(peaks[0][0],2), round(peaks[0][-1],2)
    Uc = round(p1 - p0,1)
#xnew = np.linspace(dist[0], dist[-1], num=((len(coupling)//10)+1)*10, endpoint=True)
    #arx =np.linspace(round(peaks[0][0],1),round(peaks[0][-1],1),((round(Uc,0)/10)+1,endpoint=True)
    arx =np.linspace(int(p0),int(p1),(int(Uc/10)+1),endpoint=True)
    ary = np.zeros(int(Uc/10)+1)+ylev
    ar = (arx,ary)

    #Uc = round(abs(peaks[0][0]) - abs(peaks[1][-1]),1)
    Ucstr = 'U =  ' + str(Uc) + ' mV'
    Uxy = (dips[0][0], ylev)
    U = {'U_coulomb': Uc, 'U_diagram_xy': Uxy, 'Ustring': Ucstr}

    onofn = round(ypeaks[0] / ydips[0],2) # on-on ratio neg-side
    onofp = round(ypeaks[-1] / ydips[0],2) # on-on ratio pos-side
    onofg = round(G_di['dIdV_max_global'][1] / G_di['dIdV_min_global'][1],2) # on-on ratio all
    onoff = [onofn, onofg, onofp]
    onoftxt = ['on_off_neg = '+ str(onofn), 'on_off_window = ' + str(onofg), \
               'on_off_pos = ' + str(onofp)]
    onof = (onoff, onoftxt)
    di = {'extrema': ex, 'Uline':ar, 'onoff': onof, 'U': U}
    return di

def amplitude(pd):
    #px = pd[0][0]
    py = pd[0][1]
    #dx = pd[1][0]
    dy = pd[1][1]
    if len(dy) < len(py):
        l = len(dy)
    else:
        l = len(py)
    amps = []
    for i in range(l):
        a = (py[i]-dy[i])
        amps.append(a)
    amp = np.average(amps)
    return amps, amp

def G_scale(yall, fac):
    yall[0][1],yall[1][1],yall[2][1] = \
        yall[0][1]*fac,yall[1][1]*fac,yall[2][1]*fac
    return yall

def time_limit(vlimit, varr, tarr, yarr = [0]):
    i1 = np.where(np.isclose(varr, vlimit[0], atol=6E-1))[-1]
    i2 = np.where(np.isclose(varr, vlimit[1], atol=6E-1))[-1]
    tlim = [tarr[i1][0], tarr[i2][0]]
    m = (tarr >= tlim[0]) & (tarr <= tlim[1])
    tc = tarr[m]
    vc = varr[m]
    yc = yarr[m]
    return  tlim, vc, tc, yc

def Kondo_fig_title(results, KoTemp, T_meas, Model):
    #res = results.values()
    if Model == 'slz':
        fw = results['slz_1fwhm']
        fwh = str(round(fw,3))
        am = results['slz_1amplitude']
        amp = str(round(am,1))

        ti = 'FWHM = ' + fwh + ' mV,' + ' amplitude = ' + amp + ','  \
                + ' \n ' + '${T_K} \ = \ $' + str(KoTemp) + ' K, ' + '(' + \
                    'center = ' + str(T_meas) + ' K)'
    if Model == 'gau':
        fw = results['gau_1fwhm']
        fwh = str(round(fw,3))
        am = results['gau_1amplitude']
        amp = str(round(am,1))

        ti = 'FWHM = ' + fwh + ' mV,' + ' amplitude = ' + amp + ','  \
            + ' \n ' + '${T_K} \ = \ $' + str(KoTemp) + ' K, ' + '(' + \
                'T_insitu = ' + str(T_meas) + ' K)'
    if Model == 'fano':
        fw = results['fa_1sigma']
        fwh = str(round(fw,3))
        am = results['fa_1amplitude']
        amp = str(round(am,2))
        q = results['fa_1q']
        qst = str(round(q,3))

        ti = 'FWHM = ' + fwh + ' mV,' + ' amplitude = ' + amp + ','  \
        + ' \n ' + '${T_K} \ = \ $' + str(KoTemp) + ' K, ' + '(' + \
            'T_insitu = ' + str(T_meas) + ' K),' + ' q = ' + qst
    if Model == 'fano_n-':
        fw = results['fan_1sigma']
        fwh = str(round(fw,3))
        am = results['fan_1amplitude']
        amp = str(round(am,2))
        q = results['fan_1q']
        qst = str(round(q,3))

        ti = 'FWHM = ' + fwh + ' mV,' + ' amplitude = ' + amp + ','  \
        + ' \n ' + '${T_K} \ = \ $' + str(KoTemp) + ' K, ' + '(' + \
            'T_insitu = ' + str(T_meas) + ' K),' + ' q = ' + qst
    return ti

def dct_np_tolist(odct, adkey, adlst):
    mla = {} # multi lorentzian arrays
    for i,j in enumerate(odct.values()):
        if type(j) != float and type(j) != np.float64:
            mla[list(odct.keys())[i]] = list(j)
        else:
            leng = len(list(odct.values())[i-1])
            mla[list(odct.keys())[i]] = list(np.linspace(j,j,leng))
    for k, l in enumerate(adkey):
        mla[l] = adlst[k]
    return mla

def sorter(V,y):
# =============================================================================
#     yuni = []
#     for i in range(len(y)):
#         yui = np.unique(y[i],return_index=True, return_inverse=True)[1]
#         yu = y[i][yui]
#         yuni.append(yu)
# =============================================================================
    args = np.argsort(V)
    V_ = V[args]
    yn = []
    for i in range(len(y)):
        y_ = y[i][args]
        yn.append(y_)
    return V_, yn

def ndr_calc(x,y,ra,facx=1000, facy = 1):
    m = (x >= ra[0]) & (x <= ra[1])
    xm, ym = x[m]/facx, y[m]/facy
    area = integ.simps(ym,xm,even='avg')
    ra = (ra[1]-ra[0])/facx
    area_r = area/ra
    return [xm.tolist(), ym.tolist()], [area, area_r]

def find_y_zero(x,y,rangs):
    """

    """
    r1, r2 = rangs
    y=np.array(y)
    mr1 = (x >= r1[0]) & (x <= r1[1])
    yr1 = y[mr1]

    mr2 = (x >= r2[0]) & (x <= r2[1])
    yr2 = y[mr2]
    yz1, yz2 = min(abs(yr1)), min(abs(yr2))
    v1, v2 = x[np.where(abs(y) == yz1)[0][0]], x[np.where(abs(y) == yz2)[0][0]]  #zero index 1
    return [v1, v2], [yz1, yz2]

def window(V1,y,Vlim, *V2):
    V = V1
    I,LIA1,LIA2,Iten = y[0], y[1], y[2], y[3]
    lim_ind = []
    for i,val in enumerate(V):
        if val > Vlim[0] and val < Vlim[1]:
            #print(val)
            lim_ind.append(i)
    V,I,LIA1,LIA2,Iten = V[lim_ind],I[lim_ind],LIA1[lim_ind],\
        LIA2[lim_ind],Iten[lim_ind]
    if V2:
        V2 = V2[0][lim_ind]
        return [V,V2,I,LIA1,LIA2,Iten]
    else:
        return [V,I,LIA1,LIA2,Iten]

def Vtest(V):
    Vok = []
    Vn = []
    for i in range(len(V)-1):
        if V[i]<V[i+1]:
            Vok.append((i, V[i]))
        else:
            Vn.append((i, V[i]))
            print('not in order')
    return [Vok,Vn]

def interp(x,y,fac=20):
    #xd = []
    xd = (np.linspace(min(x),max(x),len(x)*fac))
    #xd.extend(np.linspace(min(x1),max(x1),len(x1)*fac))
    yd = []
    for i in range(len(y)):
        #tck = interpolate.splrep(sorted(x), y[i], k=3)
        tck = interpolate.splrep(x, y[i], k=3)
        splev = [interpolate.splev(xd, tck, der=0)]
        yd.extend(splev)
    return [xd,yd]

def interp1d(x,y,fac=20):
    xd = (np.linspace(min(x),max(x),len(x)*fac))
    yd = []
    for i in range(len(y)):
        #tck = interpolate.splrep(sorted(x), y[i], k=3)
        tck = interpolate.interp1d(x, y[i])
        splev = [tck(xd)]
        yd.extend(splev)
    return [xd,yd]

def smooth(y,filterparams):#=[21,0,0,0.5,'interp']):
    """y can be a List of arrays"""
    window_sm = filterparams[0]
    polyorder = filterparams[1]
    deriv = filterparams[2]
    delta_sm = filterparams[3]
    mode = filterparams[4]
    Sm = []
    for i in range(len(y)):
        Smi = [signal.savgol_filter(y[i], window_sm,\
                        polyorder, deriv, delta_sm, mode=mode)]
        Sm.extend(Smi)
    return Sm

def Reflection(x,y):
    rx = [] #fpx is for flip x axis
    ry = []
    px = [] #positive x
    py = []
    nx = []
    ny = []
    for i in range(len(x)):
        if x[i] < 0 :
            rx.append(x[i]*(-1))
            ry.append(y[i]*(-1))
            nx.append(x[i])
            ny.append(y[i])
        elif x[i] > 0 :
            px.append(x[i])
            py.append(y[i])
    if len(px) != len(nx):
        if len(px)>len(nx):
            dif = len(px) - len(nx)
            px = px[:-dif]
            py = py[:-dif]
        elif len(nx) > len(px):
            dif = len(nx) - len(px)
            nx = nx[:-dif]
            ny = ny[:-dif]
            rx = rx[:-dif]
            ry = ry[:-dif]
    return [rx,ry], [px, py], [nx,ny]

def xy_shift(P, dip_n, peak_n):
    """shift spectrum in x and y direction so that Au peak is at the same
    position on pos and neg bias side. y direction is shifted so, that both
    the Au peaks have the same height. For shifting y axis the dip minimum must be
    negative.
    """
    # shift x:
    d1, p1 =P[0][1][0][dip_n], P[0][0][0][peak_n]
    shift1 = (d1+p1)/2
    if p1<d1:
        if d1 < 0:
            shift2= abs(d1 - shift1)+d1 #prooved
        else:
            if p1 >0:
                shift2 = -shift1
            else:
                shift2 = shift1
    else:
        if p1 < 0:
            shift2= abs(p1 - shift1)+p1 #prooved
        else:
            if d1 >0:
                shift2 = -shift1
            else:
                shift2 = shift1
    Vshift = shift2
    #shift y axis:
    dy1, py1 =P[0][1][1][dip_n], P[0][0][1][peak_n]
    yshift1 = (dy1+py1)/2
    if abs(py1)>abs(dy1):
        yshift2 = -abs(dy1 - yshift1)+dy1 #prooved
    else:
        yshift2 = -yshift1
    yshift = yshift2
    return [Vshift, yshift]

def derivatives(x, y):
    """y can be a list of arrays"""
    dy =[]
    for i in range(len(y)):
        #dyi = [np.gradient(x, y[i])]
        dyi = [np.gradient(y[i],x)]
        dy.extend(dyi)
    return dy

def Integrals(x, y):
    """y can be a List of arrays"""
    Y = []
    for i in range(len(y)):
        Yi = [integ.cumtrapz(y[i], x)]
        Y.extend(Yi)
    return Y

def float_to_txt(float_array,fac=1):
    P = [i*fac for i in float_array]
    Pr = [int(i) for i in P]
    txt = [str(i) for i in Pr]
    return txt #txt_array

def symmcurves(curves, Extr, SExtr, regime):
    if regime == 'IETS':
        sig = curves[0]
        E = Extr[0]
        sE = SExtr[0]
    elif regime == 'PCS':
        sig = curves[1]
        E = Extr[1]
        sE = SExtr[1]
    elif regime == 'both':
        x = np.concatenate((curves[0][0],curves[1][0]))
        y = np.concatenate((curves[0][1],curves[1][1]))
        sig = [x,y]
        Ex = np.concatenate((Extr[0][0],Extr[1][0]))
        Ey = np.concatenate((Extr[0][1],Extr[1][1]))
        E = [Ex,Ey]
        sEx = np.concatenate((SExtr[0][0],SExtr[1][0]))
        sEy = np.concatenate((SExtr[0][1],SExtr[1][1]))
        sE = [Ex,Ey]
    f = plt.figure(25, figsize = [7.68, 5.76] )
    ax = f.add_subplot(111)
    ax.scatter(sig[0],sig[1], s=1)
    ax.scatter(E[0],E[1], s=25)
    return sig, f, regime, E, sE

def Conductance(V,y,where,facy): 
    absdi = abs(V-where)
    midiff = min(absdi)
    Vi = np.where(absdi == midiff)[0]
    G = y[1][Vi]/facy
    return list(G)[0]

def smooth_param(x,sm_degree):
    """ degree of smoothin can be 'rough', 'middle', or 'fine'.
    """
    lx = len(x)
    if lx < 150:
        lx = lx + (150 - lx)
    if sm_degree == 'very_rough':
        win = round(lx/180)
        if win%2 == 0:
            win = win+1
        polyorder = 6
        deriv = 0
        delta = 0
        last = 'interp'
    if sm_degree == 'rough':
        win = round(lx/90)
        if win%2 == 0:
            win = win+1
        polyorder = 5
        deriv = 0
        delta = 0
        last = 'interp'
    elif sm_degree == 'middle':
        win = round(lx/40)
        if win%2 == 0:
            win= win+1
        polyorder = 3
        deriv = 0
        delta = 0
        last = 'interp'
    elif sm_degree == 'fine':
        win = round(lx/15)
        if win%2 == 0:
            win=win+1
        polyorder = 2
        deriv = 0
        delta = 0
        last = 'interp'
    elif sm_degree == 'very_fine':
        win = round(lx/6)
        if win%2 == 0:
            win=win+1
        polyorder = 1
        deriv = 0
        delta = 0
        last = 'interp'
    return [win, polyorder, deriv, delta, last]








### peak analysers #####################################################################
def Peak_diff(V, iets_, Pkdiff, pp, width=[1.0,1], fac = 100, lim=None): # (before peak flip)
    """yas are the as measured ones, lah --> lookahead
    """
    iets = iets_*fac
    lah, delta = pp  #findpeak parameters
    if lim != None:
        mask = (V >= lim[0]) & (V <= lim[1])
        Ex=pedet.peakdetect(iets[mask],V[mask],lookahead=lah,delta=delta) #Extrema
    else:
        Ex=pedet.peakdetect(iets,V,lookahead=lah,delta=delta)
    iets = iets/fac
    ########### unpack the list with peak positions ###########
    peaks = np.vstack((np.transpose(Ex[0])[0],np.transpose(Ex[0])[1]/fac))
    dips = np.vstack((np.transpose(Ex[1])[0],np.transpose(Ex[1])[1]/fac))
    extrema = [np.concatenate((peaks[0],dips[0])), \
        np.concatenate((peaks[1],dips[1]))]
    ###### Flip min and max extrema from negative bias region and seperate  #####
    #########            those from the positive side           ###
    fpx,fpy = Reflection(peaks[0], peaks[1])[0]  #fpx is for flip peaks x axis
    fdx, fdy = Reflection(dips[0], dips[1])[0] #for flip dips x axis
    ppx, ppy = Reflection(peaks[0], peaks[1])[1] #positive peaks x
    pdx,pdy = Reflection(dips[0], dips[1])[1]#positive dips x
    Vpos,ypos = Reflection(V,iets)[1]

    ################### Flip the smoothed curve and the Bias ##################
    [Vflip, ietsflip],[V_poshalf, iets_poshalf],[n0,n1] = Reflection(V, iets)
    ietsflipfull=[]
    Vflipfull =[]
    for i in range(len(V)):
            ietsflipfull.append(iets[i]*(-1))
            Vflipfull.append(V[i]*(-1))
### find symmetric extrema differenz between the one on neg and pos Bias side
    ###       list both the mirrored ones and the pos ones   #######
    #fpx = np.transpose(fpx)
    Sym_dip_fx = [] # symmetric dips flipped x-axis
    Sym_dip_fy = []
    Sym_peak_x = [] # symmetric peaks positive x-axis
    Sym_peak_y = []
    Sym_dip_x = []
    Sym_dip_y = []
    Sym_peak_fx = []
    Sym_peak_fy = []
    #### for peaks vs dips
    for x in range(len(fdx)):
            for i in range(len(ppx)):
                if abs(ppx[i] - fdx[x]) < Pkdiff and ppx[i] > 0.0001:
                    Sym_dip_fx.append(fdx[x])
                    Sym_dip_fy.append(fdy[x])
                    Sym_peak_x.append(ppx[i])
                    Sym_peak_y.append(ppy[i])
    #pos_x =  np.concatenate((ppx,pdx))
    #pos_y =  np.concatenate((ppy,pdy))
#====== for dips vs peaks
    for x in range(len(fpx)):
            for i in range(len(pdx)):
                if abs(pdx[i] - fpx[x]) < Pkdiff and pdx[i] > 0.0001 :
                    Sym_peak_fx.append(fpx[x])
                    Sym_peak_fy.append(fpy[x])
                    Sym_dip_x.append(pdx[i])
                    Sym_dip_y.append(pdy[i])
    #alldips_x =  np.concatenate((Sym_dip_x,Sym_peak_fx))
    #symm_alldips_y =  np.concatenate((Sym_dip_y,Sym_peak_fy))
    symm_dips_pos = [Sym_dip_x, Sym_dip_y]
    symm_peaks_pos = [Sym_peak_x,Sym_peak_y]
    #sypos_Extr = [np.array(symm_peaks_pos), np.array(symm_dips_pos)]
    sypos_Extr = [symm_peaks_pos, symm_dips_pos]
    if len(np.unique(Sym_dip_x)) < len(Sym_dip_x):
        print('/n /n peak symmetry multiple detection --> narrow diff criterium')
 ########### replot only the symmetric Peaks  ##########
    Vsympeaks  = []
    Ysympeaks = []
    Ysymdips = []
    iterable =[]
    Vsymdips = []
    Wp, Wm = width[0]*Pkdiff, width[1]*Pkdiff
    for i in range(len(Sym_peak_y)):
        iterable.clear()
        for x in range(len(V)):
            if Sym_peak_x[i]-Wm < V[x] < Sym_peak_x[i]+Wp:
                iterable.append(x)
        Vsympeaks.extend(V[iterable])
        Ysympeaks.extend(iets[iterable])
    sypecurve = np.array([Vsympeaks, Ysympeaks])

    for i in range(len(Sym_dip_y)):
        iterable.clear()
        for x in range(len(V)):
            if Sym_dip_x[i]-Wm < V[x] < Sym_dip_x[i]+Wp:
                iterable.append(x)
        Vsymdips.extend(V[iterable])
        Ysymdips.extend(iets[iterable])
    sydicurve = np.array([Vsymdips, Ysymdips])
    return [peaks,dips], extrema, sypos_Extr, [sypecurve,sydicurve]

def divid_peak_extractor(V, y, findpeak_parameters, fac):
    lah, delta = findpeak_parameters
    Extrma = pedet.peakdetect(np.array(y)*fac, V, lah, delta)
    peaks = np.transpose(Extrma[0])
    peaks_xpos, peaks_ypos = [], []
    for i in range(len(peaks[1])):
        if peaks[1][i] > 0:
            peaks_xpos.append(peaks[0][i])
            peaks_ypos.append(peaks[1][i]/fac)
    peaks = [peaks_xpos, peaks_ypos]
    dips = np.transpose(Extrma[1])
    dips_xpos, dips_ypos = [], []
    for i in range(len(dips[1])):
        if dips[1][i] < 0:
            dips_xpos.append(dips[0][i])
            dips_ypos.append(dips[1][i]/fac)
    dips = [dips_xpos, dips_ypos]
    return peaks,dips

def Peak_divide(V, iets, findpeak_parameters,lim=[0,200], sig='iets',fac=1):
    lah, delta = findpeak_parameters
    try:
        mask = (V >= lim[0]) & (V <= lim[1])
        Vm, ietsm = V[mask], iets[mask]
    except ValueError:
        Vm, ietsm = V, iets
    ################### Flip the smoothed curve and the Bias ##################
    Refl = Reflection(Vm, ietsm) #returns : Refl, pos, neg
    Vp, pos = Refl[1]
    Vn, neg = Refl[2]
    flip = Refl[0][1]
    neg.reverse()
    flip.reverse()
    #Vr = Refl[1][0]

# === new spectrum ============================0
    bkg = [(x-y)/2 for x, y in zip(pos, flip)]
    added = [(x+b)/2 for x, b in zip(pos ,flip)]
    subs = [x-b for x, b in zip(pos,bkg)]
    #add_bkg = [x-b for x, b in zip(added ,bkg)]
    peaks1, dips1 = divid_peak_extractor(Vp, subs, findpeak_parameters, fac)
    peaks2, dips2 = divid_peak_extractor(Vp, added, findpeak_parameters, fac)
    #peaks3,dips3 = divid_peak_extractor(V, add_bkg, findpeak_parameters)
    specsd = {}
    #if signal == 'didv'
    #specsd['window_spectrum_' + sig] = [V.tolist(), iets.tolist()]
    specsd['neg_spec_'+ sig] = [Vn, neg]
    specsd['pos_spec_'+ sig] = [Vp, pos]
    specsd['bkg_spec_'+ sig] = bkg # bkg -> background
    specsd['bkg_subtract_'+ sig] = subs
    specsd['antisymm_spec_'+ sig] = added
    specs = [[Vn,Vp], neg, pos,flip,bkg,subs,added]
    #peaks = np.transpose(Peaks[0])
    #dips = np.transpose(Peaks[1])
    return [peaks1,dips1],[peaks2,dips2],specs, specsd



### fit functions ########################################################################
def fit_peakfilter(Varray,extrema_x,extrema_y):
    lim_ind = []
    for i,val in enumerate(extrema_x):
        if val > min(Varray) and val < max(Varray):
            lim_ind.append(i)
    V_positions, amplitude = extrema_x[lim_ind], extrema_y[lim_ind]
    return V_positions, amplitude

def background(prefix='bkg_',x=None, shape = 'linear'):
    if shape == 'quadratic':
        model = QuadraticModel(prefix='bkg_') #f(x; a, b, c) = a x^2 + b x + c
        params = model.make_params(a=0, b=0, c=0)
    elif shape == 'linear':
        model = LinearModel(x,prefix='bkg_') #f(x; m, b) = m x + b
        params = model.make_params()
        params['bkg_slope'].set(1,min=-10,max=10)
        params['bkg_intercept'].set(0,min=-1E-6, max=15)
    elif shape=='const':
        model = ConstantModel(x, prefix='offset_') #f(x; m, b) = m x + b
        params = model.make_params()
        params['offset_c'].set(0,min=-15, max=15)
    elif shape == None:
        model = ConstantModel(x, prefix='offset_')
        params = model.make_params()
        params['offset_c'].set(0,min=-1E-9, max=1E-9)
    return model, params

def fits(x,y,model):
    params = model.guess(y, x)
    result = model.fit(y, params, x=x)
    comps = result.eval_components()
    report = result.fit_report(min_correl=0.5)
    best_fit = result.best_fit
    #gueses = (peak_positions, am)
    #print(result.fit_report(min_correl=0.5))
    return result, comps, best_fit, report, [], model, params

def pfits(x, y, extrema_x, extrema_y, bkg, prx, shape):
    prx = str(prx)
    if type(extrema_x) == np.ndarray:
        filt_ex = fit_peakfilter(x,extrema_x, extrema_y)
        peak_positions, am = filt_ex[0], filt_ex[1]
    else:
        peak_positions, am = extrema_x, extrema_y

    model, params = background(x,shape = bkg)
    for i, cen in enumerate(peak_positions):
        peak, pars = shape(prx  % (i+1), cen, am[i])
        model = model + peak
        params.update(pars)
    init = model.eval(params, x=x)
    result = model.fit(y, params, x=x)
    comps = result.eval_components()
    report = result.fit_report(min_correl=0.5)
    best_fit = result.best_fit
    gueses = (peak_positions, am)
    #print(result.fit_report(min_correl=0.5))
    return result, comps, best_fit, report, gueses, model, params

def pfits_bw(x, y, extrema_x, extrema_y, bkg, prx, shape, qlim):
    prx = str(prx)
    if type(extrema_x) == list:
        filt_ex = fit_peakfilter(x,extrema_x, extrema_y)
        peak_positions, am = filt_ex[0], filt_ex[1]
    else:
        peak_positions, am = extrema_x, extrema_y
    model, params = background(x,shape = bkg)
    for i, cen in enumerate(peak_positions):
        peak, pars = shape(prx  % (i+1), cen, am[i], qlim=qlim)
        model = model + peak
        params.update(pars)
    init = model.eval(params, x=x)
    result = model.fit(y, params, x=x)
    comps = result.eval_components()
    report = result.fit_report(min_correl=0.5)
    best_fit = result.best_fit
    gueses = (peak_positions, am)
    #print(result.fit_report(min_correl=0.5))
    return result, comps, best_fit, report, gueses, model, params

def Fit_linear(x, y):
    """other background could be quadratic or none"""
    #prefix = 'lin'
    shape= LinearModel()
    fit = fits(x, y, shape)
    return fit[0], fit[1], fit[2], fit[3], fit[4], fit[5], fit[6]

def fano_norm(x, amplitude=1.0, center=0.0, sigma=1.0, q=0):
    """Return a normalized Fano lineshape.

    breit_wigner(x, amplitude, center, sigma, q) =
        amplitude*(q*sigma/2 + x - center)**2 /
            ( (sigma/2)**2 + (x - center)**2 )

    """
    amp = amplitude
    gam = sigma/2
    return ((amp/(1+q**2))*((q*gam+x-center)**2)/((gam*gam+(x-center)**2)))

def Fano_n_Model(prefix): # ad normalized fano model
    r"""A model based on a Breit-Wigner-Fano function (see
    https://en.wikipedia.org/wiki/Fano_resonance), with four Parameters:
    ``amplitude`` (:math:`A`), ``center`` (:math:`\mu`),
    ``sigma`` (:math:`\sigma`), and ``q`` (:math:`q`).

    .. math::

        f(x; A, \mu, \sigma, q) = \frac{A/(1+q^2)* ((q\sigma/2 + x - \mu)^2}{(\sigma/2)^2 + (x - \mu)^2)}
    """
    shape = Model(fano_norm, prefix=prefix)
    #shape.set_param_hint('sigma', min=0.0)
    return shape

def add_Fano_n(prefix,center=0, amp=1E-2, sigma=10, q=0, qlim=500): # add fano model
    r"""A model based on a Breit-Wigner-Fano function (see
    https://en.wikipedia.org/wiki/Fano_resonance), with four Parameters:
    ``amplitude`` (:math:`A`), ``center`` (:math:`\mu`),
    ``sigma`` (:math:`\sigma`), and ``q`` (:math:`q`).

    .. math::

        f(x; A, \mu, \gamma, q) =
    """
    amplitude = amp
    peak = Fano_n_Model(prefix=prefix)
    pars = peak.make_params()
    pars[prefix + 'center'].set(center, min=center-1E-1, max=center+1E-1)
    pars[prefix + 'amplitude'].set(amplitude,max=10, min=1E-4)
    pars[prefix + 'sigma'].set(sigma, min=1E-3, max=100)
    pars[prefix + 'q'].set(q, min = -qlim, max = qlim)
    return peak, pars

def Fit_Fano_n(x, y, extrema_x, extrema_y, bkg , qlim):
    """other background could be quadratic or none"""
    prefix = 'fan_%d'
    shape = add_Fano_n
    fit = pfits_bw(x, y, extrema_x, extrema_y, bkg, prefix, shape, qlim=qlim)
    return fit[0], fit[1], fit[2], fit[3], fit[4], fit[5], fit[6]

def fano(x, amplitude=1.0, center=0.0, sigma=1.0, q=1.0):
    """Return a Breit-Wigner-Fano lineshape.

    breit_wigner(x, amplitude, center, sigma, q) =
        amplitude*(q*sigma/2 + x - center)**2 /
            ( (sigma/2)**2 + (x - center)**2 )

    """
    amp = amplitude
    gam = sigma/2
    return ((amp)*((q*gam+(x-center))**2)/((gam*gam+(x-center)**2)))

def Fano_Model(prefix): # ad normalized fano model
    r"""A model based on a Breit-Wigner-Fano function (see
    https://en.wikipedia.org/wiki/Fano_resonance), with four Parameters:
    ``amplitude`` (:math:`A`), ``center`` (:math:`\mu`),
    ``sigma`` (:math:`\sigma`), and ``q`` (:math:`q`).

    .. math::

        f(x; A, \mu, \sigma, q) = \frac{A/(1+q^2)* ((q\sigma/2 + x - \mu)^2}{(\sigma/2)^2 + (x - \mu)^2)}
    """
    shape = Model(fano, prefix=prefix)
    shape.set_param_hint('sigma', min=0.0)
    return shape

def add_fano(prefix,center=0, amp=1E-2, sigma=1, q=0, qlim=250):
    r"""A model based on a Breit-Wigner-Fano function (see
    https://en.wikipedia.org/wiki/Fano_resonance), with four Parameters:
    ``amplitude`` (:math:`A`), ``center`` (:math:`\mu`),
    ``sigma`` (:math:`\sigma`), and ``q`` (:math:`q`).

    .. math::

        f(x; A, \mu, \sigma, q) = \frac{A (q\sigma/2 + x - \mu)^2}{(\sigma/2)^2 + (x - \mu)^2}

    """
    amplitude = amp
    peak = Fano_Model(prefix=prefix)
    pars = peak.make_params()
    pars[prefix + 'center'].set(center, min=center-1E-1, max=center+1E-1)
    pars[prefix + 'amplitude'].set(amplitude,\
                max=10, min=1E-4)
    pars[prefix + 'sigma'].set(sigma, min=1E-3, max=100)
    pars[prefix + 'q'].set(q, min = -qlim, max = qlim)
    return peak, pars

def Fit_Fano(x, y, extrema_x, extrema_y, bkg = 'linear',qlim=10):
    """other background could be quadratic or none"""
    prefix = 'fa_%d'
    shape = add_fano
    fit = pfits_bw(x, y, extrema_x, extrema_y, bkg, prefix, shape, qlim)
    return fit[0], fit[1], fit[2], fit[3], fit[4], fit[5], fit[6]


def add_peak_BW(prefix,center=0, amp=1E-2, sigma=1, q=0, qlim=250):
    peak = BreitWignerModel(prefix=prefix)
    pars = peak.make_params()
    pars[prefix + 'center'].set(center, min=center-1E-1, max=center+1E-1)
    pars[prefix + 'amplitude'].set(amplitude,\
                max=10, min=1E-4)
    pars[prefix + 'sigma'].set(sigma, min=1E-3, max=100)
    pars[prefix + 'q'].set(q, min = -qlim, max = qlim)   
    return peak, pars

def Fit_BreitWigner(x, y, extrema_x, extrema_y, bkg = 'linear',qlim=2):
    """other background could be quadratic or none"""
    prefix = 'bw_%d'
    shape = add_peak_BW
    fit = pfits_bw(x, y, extrema_x, extrema_y, bkg, prefix, shape, qlim)
    return fit[0], fit[1], fit[2], fit[3], fit[4], fit[5], fit[6]

def add_peak_lor(prefix, center=1, height=1, sigma=20):
    peak = LorentzianModel(prefix=prefix)
    pars = peak.make_params()
    pars[prefix + 'center'].set(center, min=center-7, max=center+7)
    pars[prefix + 'height'].set(height,\
                  max=height+(abs(height)/5),\
                      min=height-(abs(height)/5))
    pars[prefix + 'sigma'].set(sigma, min=1E-1, max=250)
    return peak, pars

def Fit_lorentz_si(x, y, extrema_x, extrema_y, bkg = 'linear'):
    """other background could be quadratic, const or none"""
    prefix = 'lz_%d'
    shape= add_peak_lor
    fit = pfits(x, y, extrema_x, extrema_y, bkg, prefix, shape)
    return fit[0], fit[1], fit[2], fit[3], fit[4], fit[5], fit[6]

def add_peak_Gau(prefix, center=1, height=1, sigma=20):
    peak = GaussianModel(prefix=prefix)
    pars = peak.make_params()
    pars[prefix + 'center'].set(center, min=center-5, max=center+5)
    pars[prefix + 'height'].set(height,\
                  max=height+(abs(height)/9),\
                      min=height-(abs(height)/9))
    pars[prefix + 'sigma'].set(sigma, min=1E-3, max=200)
    return peak, pars

def Fit_Gauss(x, y, extrema_x, extrema_y, bkg = 'linear'):
    """other background could be quadratic or none"""
    prefix = 'gau_%d'
    shape= add_peak_Gau
    fit = pfits(x, y, extrema_x, extrema_y, bkg, prefix, shape)
    return fit

def add_peak_slor(prefix, center=1, height=1, sigma=15, sigma_r=25):
  peak = SplitLorentzianModel(prefix=prefix)
  pars = peak.make_params()
  pars[prefix + 'center'].set(center, min=center-15, max=center+15)
  pars[prefix + 'height'].set(height,\
                  max=height+(abs(height)/10),\
                      min=height-(abs(height)/10))
  pars[prefix + 'sigma'].set(sigma, min=1E-2, max=1E2)
  pars[prefix + 'sigma_r'].set(sigma, min=1E-2, max=1E2)
  return peak, pars

def Fit_slor(x, y, extrema_x, extrema_y, bkg = 'linear'):
    """other background could be quadratic or none"""
    prefix = 'slz_%d'
    shape= add_peak_slor
    fit = pfits(x, y, extrema_x, extrema_y, bkg, prefix, shape)
    return fit[0], fit[1], fit[2], fit[3], fit[4], fit[5], fit[6]

def add_peak_sgau(prefix,center=1,amplitude=1,sigma=15,gamma=20):
    peak = SkewedGaussianModel(prefix=prefix)
    pars = peak.make_params()
    pars[prefix + 'center'].set(center, min=center-1, max=center+1)
    pars[prefix + 'amplitude'].set(amplitude, max=amplitude+\
                       (abs(amplitude)/10),min=amplitude-(abs(amplitude)/10))
    pars[prefix + 'sigma'].set(sigma, min=1E-5, max=70)
    pars[prefix + 'gamma'].set(gamma, min=1E-5, max=70)
    #pars[prefix + 'skew'].set(skew, min=0, max=1E-3)
    return peak, pars

def Fit_sGauss(x, y, extrema_x, extrema_y, bkg = 'linear'):
    """other background could be quadratic or none"""
    prefix = 'sgau_%d'
    shape= add_peak_sgau
    extrema_y_ = []
    for i,val in enumerate(extrema_y):
        new = val/np.sqrt(2*np.pi)
        extrema_y_.append(new)
    exy = np.array(extrema_y_)
    fit = pfits(x, y, extrema_x, extrema_y, bkg, prefix, shape)
    return fit[0], fit[1], fit[2], fit[3], fit[4], fit[5], fit[6]

def add_peak_voi(prefix, center=1, height=1, sigma=15, gamma=20):
    peak = VoigtModel(prefix=prefix)
    pars = peak.make_params()
    pars[prefix + 'center'].set(center, min=center-5, max=center+5)
    pars[prefix + 'height'].set(height, max=height+(abs(height)/10),\
                        min=height-(abs(height)/10))
    pars[prefix + 'sigma'].set(sigma, min=1E-5, max=75)
    pars[prefix + 'gamma'].set(gamma, min=1E-5, max=75)
    return peak, pars

def Fit_voigt(x, y, extrema_x, extrema_y, bkg = 'linear'):
    """other background could be quadratic or none"""
    prefix = 'voi_%d'
    shape= add_peak_voi
    fit = pfits(x, y, extrema_x, extrema_y, bkg, prefix, shape)
    return fit[0], fit[1], fit[2], fit[3], fit[4], fit[5], fit[6]

def add_peak_svoi(prefix,center=0.01,amplitude=10,sigma=10,gamma=20,skew=1):
    r"""A variation of the Skewed Gaussian, this applies the same skewing
      to a Voigt distribution (see
      https://en.wikipedia.org/wiki/Voigt_distribution).  It has Parameters
      ``amplitude`` (:math:`A`), ``center`` (:math:`\mu`), ``sigma``
      (:math:`\sigma`), and ``gamma`` (:math:`\gamma`), as usual for a
      Voigt distribution, and add a Parameter ``skew``.
    .. math::
       f(x; A, \mu, \sigma, \gamma, \rm{skew}) = {\rm{Voigt}}(x; A, \mu, \sigma, \gamma)
       \Bigl\{ 1 +  {\operatorname{erf}}\bigl[
       \frac{{\rm{skew}}(x-\mu)}{\sigma\sqrt{2}}
       \bigr] \Bigr\}
    where :func:`erf` is the error function.
    """
    peak = SkewedVoigtModel(prefix=prefix)
    pars = peak.make_params()
    pars[prefix + 'center'].set(center, min=center-5, max=center+5)
    pars[prefix + 'amplitude'].set(amplitude, max=amplitude+\
                       (abs(amplitude)/10),min=amplitude-(abs(amplitude)/10))
    pars[prefix + 'sigma'].set(sigma, min=1E-5, max=70)
    pars[prefix + 'gamma'].set(gamma, min=1E-5, max=100)
    pars[prefix + 'skew'].set(skew, min=1E-7, max=50)
    return peak, pars

def Fit_svoigt(x, y, extrema_x, extrema_y, bkg = 'linear'):
    """other background could be quadratic or none"""
    prefix = 'svoi_%d'
    shape = add_peak_svoi
    fit = pfits(x, y, extrema_x, extrema_y, bkg, prefix, shape)
    return fit[0], fit[1], fit[2], fit[3], fit[4], fit[5], fit[6]

def add_peak_DHO(prefix, center=0.01, height=10, sigma=0.01, gamma=0.01):
    r"""A model based on a variation of the Damped Harmonic Oscillator (see
    https://en.wikipedia.org/wiki/Harmonic_oscillator), following the
    definition given in DAVE/PAN (see https://www.ncnr.nist.gov/dave/) with
    four Parameters: ``amplitude`` (:math:`A`), ``center`` (:math:`\mu`),
    ``sigma`` (:math:`\sigma`), and ``gamma`` (:math:`\gamma`).
    In addition, parameters ``fwhm`` and ``height`` are included as constraints
    to report estimates for full width at half maximum and maximum peak height,
    respectively.

    .. math::

        f(x; A, \mu, \sigma, \gamma) = \frac{A\sigma}{\pi [1 - \exp(-x/\gamma)]}
                \Big[ \frac{1}{(x-\mu)^2 + \sigma^2} - \frac{1}{(x+\mu)^2 + \sigma^2} \Big]

    where :math:`\gamma=kT` k is the Boltzmann constant in :math:`evK^-1`
    and T is the temperature in :math:`K`.

    """
    peak = DampedHarmonicOscillatorModel(prefix=prefix)
    pars = peak.make_params()
    pars[prefix + 'center'].set(center, min=center-5, max=center+5)
    pars[prefix + 'height'].set(height,\
                    max=height+(abs(height)/10),\
                        min=height-(abs(height)/10))
    pars[prefix + 'sigma'].set(sigma, min=1E-5, max=75)
    pars[prefix + 'gamma'].set(gamma, min=1E-5, max=75)
    return peak, pars

def Fit_Damped_Osc(x, y,extrema_x, extrema_y, bkg = 'linear'):
    """other background could be quadratic or none"""
    prefix = 'dho_%d'
    shape= add_peak_DHO
    fit = pfits(x, y, extrema_x, extrema_y, bkg, prefix, shape)
    return fit[0], fit[1], fit[2], fit[3], fit[4], fit[5], fit[6]

def add_peak_BEinst(prefix, center=0.01, amplitude=5, kt=0.01):
    """    Comments:

    - ``kt`` should be defined in the same units as ``x`` (:math:`k_B = 8.617\times10^{-5}` eV/K).
    - set :math:`kt<0` to implement the energy loss convention common in
      scattering research.
      """
    peak = ThermalDistributionModel(prefix=prefix, form='bose')
    pars = peak.make_params()
    pars[prefix + 'center'].set(center, min=center-5, max=center+5)
    pars[prefix + 'amplitude'].set(amplitude,\
                    max=amplitude+(abs(amplitude)/10),\
                        min=amplitude-(abs(amplitude)/10))
    pars[prefix + 'kt'].set(kt, min=1E-4, max=75)
    return peak, pars

def Fit_BEinst(x, y,extrema_x, extrema_y, bkg = 'linear'):
    """other background could be quadratic or none"""
    prefix = 'BE_%d'
    shape= add_peak_BEinst
    fit = pfits(x, y, extrema_x, extrema_y, bkg, prefix, shape)
    return fit[0], fit[1], fit[2], fit[3], fit[4], fit[5], fit[6]

def add_peak_Doniach(prefix, center=0.01, height=1, sigma=0.05, gamma=0.01):
    r"""A model of an Doniach Sunjic asymmetric lineshape
    (see https://www.casaxps.com/help_manual/line_shapes.htm), used in
    photo-emission, with four Parameters ``amplitude`` (:math:`A`),
    ``center`` (:math:`\mu`), ``sigma`` (:math:`\sigma`), and ``gamma``
    (:math:`\gamma`).
    In addition, parameter ``height`` is included as a constraint.

    .. math::

        f(x; A, \mu, \sigma, \gamma) = \frac{A}{\sigma^{1-\gamma}}
        \frac{\cos\bigl[\pi\gamma/2 + (1-\gamma)
        \arctan{((x - \mu)}/\sigma)\bigr]} {\bigr[1 + (x-\mu)/\sigma\bigl]^{(1-\gamma)/2}}

    """
    peak = DoniachModel(prefix=prefix)
    pars = peak.make_params()
    pars[prefix + 'center'].set(center, min=center-5, max=center+5)
    pars[prefix + 'height'].set(height, max=height+(abs(height)/10),\
                        min=height-(abs(height)/10))
    pars[prefix + 'sigma'].set(sigma, min=1E-5, max=75)
    pars[prefix + 'gamma'].set(gamma, min=1E-5, max=75)
    return peak, pars

def Fit_Doniach(x, y, extrema_x, extrema_y, bkg = 'linear'):
    """other background could be const, quadratic or none"""
    prefix = 'Do_%d'
    shape= add_peak_Doniach
    fit = pfits(x, y, extrema_x, extrema_y, bkg, prefix, shape)
    return fit[0], fit[1], fit[2], fit[3], fit[4], fit[5], fit[6]

def add_peak_logn(prefix, center=0.01, height=5, sigma=0.02):
    peak = LognormalModel(prefix=prefix)
    pars = peak.make_params()
    pars[prefix + 'center'].set(center, min=center-5, max=center+5)
    pars[prefix + 'height'].set(height,\
                  max=height+(abs(height/10)),\
                      min=height-(abs(height)/10))
    pars[prefix + 'sigma'].set(sigma, min=1E-5, max=75)
    return peak, pars

def Fit_logNormal(x, y, extrema_x, extrema_y, bkg = 'linear'):
    """other background could be quadratic or none"""
    prefix = 'logN_%d'
    shape= add_peak_logn
    fit = pfits(x, y, extrema_x, extrema_y, bkg, prefix, shape)
    return fit[0], fit[1], fit[2], fit[3], fit[4], fit[5], fit[6]




#### Figures ###########################################################################



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

def fig_ad_sub(V,y,fig,ax,plotpos,x_label,y_label,cl=[]):
    """ads a subplot to figure 'fif' at postition 'plotpos'
    plotpos is given as triple 'eg. 111' or [111] with the axis name ax
    cl:label of the curve if not then nothin happens
    """
    ax = fig.add_subplot(plotpos)
    c = ['k','b','m','c']
    a = [0.8,1,1,1]
    for i in range(len(V)):
        if len(cl) > 0:
            ax.plot(V[i],y[i],color=c[i], alpha=a[i], label=cl[i])
            handles, labels = ax.get_legend_handles_labels()
            plt.legend()
        else:
            ax.plot(V[i],y[i],color=c[i], alpha=a[i])
    ax.set_xlabel(x_label, fontsize = 16)
    ax.set_ylabel(y_label, fontsize = 16)
    ax.tick_params(axis='both' , labelsize = 12.0)
    return ax

def fig1(V,y,temp, sign, ylab):
    """ signal 'iets' or 'dG/dV' """
    #V = V*1000
    plt.close('all')
    I,LIA1,LIA2, Iten = y[0], y[1], y[2], y[3]
    f1 = plt.figure(1,figsize=(6.5,8.5))
    ax3 = f1.add_subplot(313)
    ax1 = f1.add_subplot(311, sharex=ax3)
    ax2 = f1.add_subplot(312, sharex=ax3)
    #ax4 = ax3.twinx()
    f1.tight_layout()
    ax1.plot(V,I, c='k')
    ax1.set_ylabel('current [\u00b5A]', fontsize = 16)
    ax1.tick_params(axis='both' , labelsize = 12.0)
    ax2.plot(V, LIA1, c='k')
    ax2.set_ylabel(ylab, fontsize = 16)
    ax2.tick_params(axis='both' , labelsize = 12.0)
    if sign == 'dgdv':
        ax3.plot(V, LIA2, c='k')#, label='dG/dV [G0]')
        ax3.set_ylabel('$dG/dV \ [(G/G_{0})*V^{-1}]$', fontsize = 16)
    elif sign == 'iets':
        ax3.plot(V, Iten, c='k', alpha=1)
        ax3.set_ylabel('$Intensity \ [V^{-1}]$', fontsize = 16)
    ax3.set_xlabel('Bias [mV]', fontsize = 16)
    ax3.tick_params(axis='both' , labelsize = 12.0)
    ax3.legend(['T = ' + str(temp)+ ' K'])
    f1.tight_layout()
    return f1

def fig1_all(V,y,temp, ylab):
    """ signal 'iets' or 'dG/dV' """
    #V = V*1000
    plt.close('all')
    I,LIA1,LIA2, Iten = y[0], y[1], y[2], y[3]
    f1 = plt.figure(1,figsize=(5.6,9.8))
    ax4 = f1.add_subplot(414)
    ax3 = f1.add_subplot(413, sharex=ax4)
    ax1 = f1.add_subplot(411, sharex=ax4)
    ax2 = f1.add_subplot(412, sharex=ax4)
    #ax4 = ax3.twinx()
    f1.tight_layout()
    ax1.plot(V,I, c='k')
    ax1.set_ylabel('current [\u00b5A]', fontsize = 16)
    ax1.tick_params(axis='both' , labelsize = 12.0)
    ax2.plot(V, LIA1, c='k')
    ax2.set_ylabel(ylab, fontsize = 16)
    ax2.tick_params(axis='both' , labelsize = 12.0)
    ax3.plot(V, LIA2, c='k')#, label='dG/dV [G0]')
    ax3.set_ylabel('$dG/dV \ [(G/G_{0})*V^{-1}]$', fontsize = 16)
    ax3.tick_params(axis='both' , labelsize = 12.0)
    ax4.plot(V, Iten, c='k', alpha=1)
    ax4.set_ylabel('$Intensity \ [V^{-1}]$', fontsize = 16)
    ax4.set_xlabel('Bias [mV]', fontsize = 16)
    ax4.tick_params(axis='both' , labelsize = 12.0)
    ax4.legend(['T = ' + str(temp)+ ' K'])
    f1.tight_layout()
    return f1

def Peak_fig_title(results, Modl):
    #res = results.values()
    if Modl == 'slz':
        fw = results['slz_1fwhm']
        fwh = str(round(fw,3))
        am = results['slz_1amplitude']
        amp = str(round(am,1))
        ce = results['slz_1center']
        cen = str(round(ce,1))
        of = results['bkg_intercept']
        ofs = str(round(of,1))
        ti = 'FWHM = ' + fwh + ' mV,' + ' amplitude = ' + amp + ','  \
                + ' \n ' + ' center = ' + cen + ' mV, ' + '(' + \
                    'offset = ' + ofs + ' G)'
    if Modl == 'gau' or Modl == 'lz':
        if Modl == 'gau':
            fw = results['gau_1fwhm']
            am = results['gau_1amplitude']
            ce = results['gau_1center']
            of = results['bkg_intercept']
        if  Modl == 'lz':
            fw = results['lz_1fwhm']
            am = results['lz_1amplitude']
            ce = results['lz_1center']
            of = results['bkg_intercept']
        fwh = str(round(fw,3))
        amp = str(round(am,1))
        cen = str(round(ce,1))
        ofs = str(round(of,1))
        ti = 'FWHM = ' + fwh + ' mV,' + ' amplitude = ' + amp + ','  \
                + ' \n ' + ' center = ' + cen + ' mV, ' + '(' + \
                        'offset = ' + ofs + ' G)'
    return ti

def fig_up_down(V_all, y_all, plots,posneg,cl,gylab):
    """ plots is array with the number of plots and what should be plotted.
    it can be ['I','dI/dV', 'dG/dV'] (iets or dG/dv) or just one or two of
    them. The last entry in the list is the lower plot. cl is the label
    of the curves if there are more than one in a plot.
    """
    y,yn,yp = y_all
    I,LIA1,LIA2, iets = y[0], y[1], y[2], y[3]
    In,LIA1n,LIA2n,ietsn = yn[0], yn[1], yn[2], yn[3]
    Ip,LIA1p,LIA2p,ietsp = yp[0], yp[1], yp[2], yp[3]
    V, Vn, Vp = V_all
    if len(posneg) ==2:
        Is = [I,In,Ip]
        LIA1s = [LIA1, LIA1n, LIA1p]
        LIA2s = [LIA2, LIA2n, LIA2p]
        ietss = [iets, ietsn, ietsp]
        Vs = [V, Vn, Vp]
    elif len(posneg) ==1:
        if posneg[0] == 'pos':
            Is = [I,Ip]
            LIA1s = [LIA1, LIA1p]
            LIA2s = [LIA2, LIA2p]
            ietss = [iets, ietsp]
            Vs = [V, Vp]
        elif posneg[0] == 'neg':
            Is = [I,In]
            LIA1s = [LIA1, LIA1n]
            LIA2s = [LIA2, LIA2n]
            ietss = [iets, ietsn]
            Vs = [V, Vn]
    if len(plots)==3:
        f_ud = plt.figure(1,figsize=(6,8.5))
        fig_ad_sub(Vs,Is,f_ud,'ax1',311,'Bias [mV]','current [\u00b5A]',cl=cl)
        fig_ad_sub(Vs,LIA1s,f_ud,'ax2',312,'Bias [mV]', gylab)
        if plots[2]== 'iets':
            fig_ad_sub(Vs,ietss,f_ud,'ax3',313,'Bias [mV]','IETS [1/V]')
        elif plots[2] == 'dgdv':
            fig_ad_sub(Vs,LIA2s,f_ud,'ax3',313,'Bias [mV]', \
                       '$dG/dV \ [(G/G_{0})*V^{-1}]$')
    elif len(plots)==1:
        f_ud = plt.figure(1,figsize=[6.4, 4.8])
        if plots[0] == 'I':
            fig_ad_sub(Vs,Is,f_ud,'ax1',111,'Bias [mV]','current [\u00b5A]')
        elif plots[0] == 'didv':
            fig_ad_sub(Vs,LIA1s,f_ud,'ax1',111,'Bias [mV]', gylab)
        elif plots[0]== 'iets':
            fig_ad_sub(Vs,ietss,f_ud,'ax1',111,'Bias [mV]','iets [1/V]')
        elif plots[0] == 'dgdv':
            fig_ad_sub(Vs,LIA2s,f_ud,'ax1',111,\
                       'Bias [mV]','$dG/dV [(G/G_{0})*V^{-1}]$')
    elif len(plots)==4:
        f_ud = plt.figure(1,figsize=[5.6, 9.8])
        fig_ad_sub(Vs,Is,f_ud,'ax1',411,'Bias [mV]','current [\u00b5A]',cl=cl)
        fig_ad_sub(Vs,LIA1s,f_ud,'ax2',412,'Bias [mV]', gylab)
        fig_ad_sub(Vs,LIA2s,f_ud,'ax1',413,\
                       'Bias [mV]','$dG/dV [(G/G_{0})*V^{-1}]$')
        fig_ad_sub(Vs,ietss,f_ud,'ax3',414,'Bias [mV]','IETS [1/V]')
    f_ud.tight_layout()
    return f_ud

def fig_scatter(V_all,y_all,plots,posneg):
    y,yn,yp = y_all
    I,LIA1,LIA2, iets = y[0], y[1], y[2], y[3]
    In,LIA1n,LIA2n,ietsn = yn[0], yn[1], yn[2], yn[3]
    Ip,LIA1p,LIA2p,ietsp = yp[0], yp[1], yp[2], yp[3]
    V, Vn, Vp = V_all
    f_sc = plt.figure(1,figsize=[6.4, 4.8])
    if len(posneg) ==2:
        # 's' in the following as the plural:
        Is = [I,In,Ip]
        LIA1s = [LIA1, LIA1n, LIA1p]
        LIA2s = [LIA2, LIA2n, LIA2p]
        ietss = [iets, ietsn, ietsp]
        Vs = [V, Vn, Vp]
    elif len(posneg) ==1:
        if posneg[0] == 'pos':
            Is = [I,Ip]
            LIA1s = [LIA1, LIA1p]
            LIA2s = [LIA2, LIA2p]
            ietss = [iets, ietsp]
            Vs = [V, Vp]
        elif posneg[0] == 'neg':
            Is = [I,In]
            LIA1s = [LIA1, LIA1n]
            LIA2s = [LIA2, LIA2n]
            ietss = [iets, ietsn]
            Vs = [V, Vn]
    if plots[0] == 'I':
        ad_sub_scatter(Vs,Is,f_sc,'ax1',111,'Bias [mV]','current [\u00b5A]')
    elif plots[0] == 'didv':
        ad_sub_scatter(Vs,LIA1s,f_sc,'ax1',111,'Bias [mV]', 'dI/dV [G0]')
    elif plots[0]== 'iets':
        ad_sub_scatter(Vs,ietss,f_sc,'ax1',111,'Bias [mV]','$IETS \ [V^{-1}]$')
    elif plots[0] == 'dgdv':
        ad_sub_scatter(Vs,LIA2s,f_sc,'ax1',111,'Bias [mV]',\
                       '$dG/dV \ [(G/G0)*V^{1}]$')
    elif plots[0] == 'mdidv':
        ad_sub_scatter(Vs,LIA1s,f_sc,'ax1',111,'Bias [mV]','dI/dV [mG/G0]')
    f_sc.tight_layout()
    return f_sc

def fig_compare(V,V2,y,y2,labels,fn,gyl,fp=None, dgdv = 'No',c='tab:orange'):
    """fp -> plot filter parameters legend, fn = fignumber"""
    I ,LIA1, LIA2, Iten = y
    I_2 ,LIA1_2, LIA2_2, Iten_2 = y2

    f2 = plt.figure(2, figsize=(6,8.5))
    ax3 = f2.add_subplot(313)
    ax1 = f2.add_subplot(311, sharex=ax3)
    ax2 = f2.add_subplot(312, sharex=ax3)

    ax1.plot(V,I, '--', c=c)
    ax2.plot(V, LIA1, '--', c=c)
    if fp != None:
        ax1.plot(V2,I_2,label='parameters: ' +fp)
        ax1.legend()
    else:
        ax1.plot(V2,I_2,c='k')
    ax2.plot(V2,LIA1_2, c='k')
    if dgdv == 'yes':
        #ax4 = ax3.twinx()
        ax3.plot(V,LIA2, '--', c=c )
        ax3.plot(V2, LIA2_2,'--', c='k', lw=1,label=labels[1])
        ax3.set_ylabel('dG/dV [G0]', fontsize = 16)
    elif dgdv == 'No':
        ax3.plot(V, Iten, '--', c=c, label=labels[0])
        ax3.plot(V2, Iten_2, c='k', label=labels[1])
        ax3.set_ylabel('IETS Intensity [1/V]', fontsize = 16)
    ax1.set_ylabel('current [\u00b5A]', fontsize = 16)
    ax1.tick_params(axis='both' , labelsize = 12.0)
    ax2.set_ylabel(gyl, fontsize = 16)
    ax2.tick_params(axis='both' , labelsize = 12.0)
    ax3.set_xlabel('Bias [mV]', fontsize = 16)
    ax3.tick_params(axis='both' , labelsize = 12.0)
    ax3.legend()
    f2.tight_layout()
    return f2

def peakPlot(x, y, pd, spec, ylab, meth='ansy', show='both', xlab='Bias [mV]'):
    """ pd are peaks and dips in form of a list:
        --> all peaks and dips: [peaks,dips](Peak_flip[0]),
        symmetric peaks and dips: [Symmpeak,symmdip] (Peak_flip[2])
        spectum is either symmetry or whole
        spec (spectrum): either 'whole' or 'symmetry'
        meth --> method ansy or bkg
        """
    #plt.close('all')
    if ylab == 'iets':
        ylab = '$Intensity [V^{-1}]$'
    elif ylab == 'dgdv':
        ylab = '$dG/dV \  [(G/G_{0})V^{-1}]$'
    elif ylab == 'didv':
        ylab = '$dI/dV \ [G/G_{0}]$'
    elif ylab == 'mdidv':
        ylab = '$dI/dV \ [mG/G_{0}]$'
    else:
        ylab = ylab
    R = Reflection(x,y)
    if spec == "whole":
        peaks, dips = pd
        ptxt, dtxt = float_to_txt(peaks[0]), float_to_txt(dips[0])
        fig = plt.figure(25, figsize = [7.68, 5.76] )
        ax11 = fig.add_subplot(111)
        ax11.plot(peaks[0],peaks[1], 'vr')
        ax11.plot(dips[0], dips[1] , 'vg')
        #ax11.plot(x, y, ':' , c='c' , linewidth= 2.2)
        ax11.plot(x, y, c='k' , linewidth= 1.2)
        for label, x, y in zip(ptxt, peaks[0],peaks[1]):
            ax11.annotate(label, xy = (x , y), xytext=(10, 5),color='r',\
            textcoords = 'offset points', ha='right', va='bottom',\
            fontsize = 14)
        for label, x, y in zip(dtxt, dips[0], dips[1]):
            ax11.annotate(label, xy = (x , y), xytext=(10, 5), color='g',\
            textcoords = 'offset points', ha='right', va='bottom',\
            fontsize = 14)
        ax11.set_xlabel(xlab,fontsize=18)
        ax11.set_ylabel(ylab,fontsize=18)
        ax11.tick_params(axis='both' , labelsize = 16.0)
        fig.tight_layout()
    elif spec == 'symmetry':
        peaks, dips = pd[0]
        ptxt, dtxt = float_to_txt(peaks[0]), float_to_txt(dips[0])
        specs = pd[3]
        pos = R[1]
        flip = R[0]       
        fig = plt.figure(26, figsize = [7.68 , 5.76] )
        ax21 = fig.add_subplot(111)
        ax21.plot(pos[0], pos[1], c='k' , linewidth= 1.5 , label='as measured') # ,  ,
        ax21.plot(flip[0] , flip[1] , '--r' , lw=0.8 ,\
                  label='reflected')
        if ylab == '$Intensity [V^{-1}]$' :
            ansy_spec = specs['antisymm_spec_iets']
            bkg = specs['bkg_spec_iets']
            bkgs = specs['bkg_subtract_iets']
            spx = specs['pos_spec_iets'][0]
        elif ylab == '$dG/dV \  [(G/G_{0})V^{-1}]$': 
            ansy_spec = specs['antisymm_spec_dgdv']
            bkg = specs['bkg_spec_dgdv']
            bkgs = specs['bkg_subtract_dgdv']
            spx = specs['pos_spec_dgdv'][0]
        if meth == 'ansy':
            ax21.plot(spx, ansy_spec, c='b', lw=1.5, label='antisymm')                
        elif meth == 'bkg':
            ax21.plot(spx,bkg,':', c='gray',lw=0.9,label='bkgnd')
            ax21.plot(spx,bkgs, c='b', lw=1.5, label='bkgnd_subtracted')
            
        if show == 'both':
            ax21.plot(peaks[0] , peaks[1], 'vr', label='IETS peaks')
            ax21.plot(dips[0] , dips[1], 'vg', label='PCS dips')
            for label, x, y in zip(ptxt, peaks[0] , peaks[1]):
                ax21.annotate(label, xy = (x , y), xytext=(10, 5),
                textcoords = 'offset points', ha='right', va='bottom', color='r' ,
                fontsize = 14)
            for label, x, y in zip(dtxt, dips[0] , dips[1]):
                ax21.annotate(label, xy = (x , y), xytext=(10, 5),
                textcoords = 'offset points', ha='right', va='bottom', color='g' ,
                fontsize = 14)
        elif show == 'iets':
            ax21.plot(peaks[0] , peaks[1], 'vr', label='IETS peaks')
            for label, x, y in zip(ptxt, peaks[0] , peaks[1]):
                ax21.annotate(label, xy = (x , y), xytext=(10, 5),
                textcoords = 'offset points', ha='right', va='bottom', color='r' ,
                fontsize = 14)
        elif show == 'pcs':
            ax21.plot(dips[0] , dips[1], 'vg', label='PCS dips')
            for label, x, y in zip(dtxt, dips[0] , dips[1]):
                ax21.annotate(label, xy = (x , y), xytext=(10, 5),
                textcoords = 'offset points', ha='right', va='bottom', color='g' ,
                fontsize = 14)
        ax21.set_xlabel(xlab,fontsize=18)
        ax21.set_ylabel(ylab,fontsize=18)
        ax21.tick_params(axis='both' , labelsize = 16.0)
        fig.tight_layout()
        ax21.legend()
    return fig, R[1]

# =============================================================================
# def delta_peak_plot(x,y,pd,pdd, ylab,xlabel,legend,lp):
#     """ lp is the position of the legend (lower_left, upper right,..)
#     """
#     peaks, dips = pd
#     pdif, ddif = pdd
#     pdpos,ddpos = [], []
#     for i in range(len(pdif)):
#         pdpos.append((pdif[i]+peaks[0][i]))
#     for i in range(len(ddif)):
#         ddpos.append((ddif[i]+dips[0][i]))
#     pdtxt = [str(int(round(i))) for i in pdif]
#     ddtxt = [str(int(round(i))) for i in ddif]
#     ptxt = [str(int(round(i))) for i in peaks[0]]
#     dtxt = [str(int(round(i))) for i in dips[0]]
#     fig = plt.figure()#, figsize = [7.68, 5.76] )
#     ax11 = fig.add_subplot(111)
#     ax11.plot(peaks[0],peaks[1], 'vr')
#     ax11.plot(dips[0], dips[1] , 'vg')
#     #ax11.plot(x, y, ':' , c='c' , linewidth= 2.2)
#     ax11.plot(x, y, c='k' , linewidth= 1.1)
#     ax11.scatter(x, y, c='b' , s= 2)
#     fig.tight_layout()
#     for label, x, y in zip(ptxt, peaks[0],peaks[1]):
#         ax11.annotate(label, xy = (x , y), xytext=(10, 5),color='r',\
#         textcoords = 'offset points', ha='right', va='bottom',\
#         fontsize = 14)
#     for label, x, y in zip(dtxt, dips[0], dips[1]):
#         ax11.annotate(label, xy = (x , y), xytext=(10, 5), color='g',\
#         textcoords = 'offset points', ha='right', va='bottom',\
#         fontsize = 14)
#     fig.canvas.draw()
#     fig.tight_layout()
#     for label, x, y in zip(pdtxt, pdpos, peaks[1]):
#         ax11.annotate(label, xy = (x , y), xytext=(-15, 20),color='b',\
#         textcoords = 'offset points', ha='right', va='bottom',\
#         fontsize = 15)
#     fig.tight_layout()
#     for label, x, y in zip(ddtxt, ddpos, dips[1]):
#         ax11.annotate(label, xy = (x , y), xytext=(-15, -10), color='m',\
#         textcoords = 'offset points', ha='right', va='bottom',\
#         fontsize = 15)
#     fig.canvas.draw()
#     fig.tight_layout()
#     ax11.set_xlabel(xlabel,fontsize=18)
#     ax11.set_ylabel(ylab,fontsize=18)
#     ax11.legend([legend],loc= lp)
#     ax11.tick_params(axis='both' , labelsize = 15.0)
#     fig.canvas.draw()
#     fig.tight_layout()
#     return fig
# =============================================================================

def fig_deriv(V, y, dy, y_int, legend = 'yes', dgdv='No', lim = 1):
    """if an integration is calculated"""
    I ,LIA1, LIA2, Iten = y
    dI, dLIA1, dLIA2, dIten = dy
    I_int, LIA1_int, LIA2_int, iten_int = y_int
    plt.close('all')
    f2 = plt.figure(2, figsize=(6,8.5))
    ax3 = f2.add_subplot(313)
    ax1 = f2.add_subplot(311, sharex=ax3)
    ax2 = f2.add_subplot(312, sharex=ax1)
    ax11 = ax1.twinx()
    ax22 = ax2.twinx()
    ax31 = ax3.twinx()
    ax1.plot(V, I, c='k')
    ax2.plot(V, LIA1, c='k')
    ax11.plot(V[:-1], LIA1_int, '--', c='c',  label='integrated')
    if dgdv == 'yes':
        ax22.plot(V[:-1],LIA2_int, '--', c='c', label='integrated')
        ax3.plot(V, LIA2, c='k')
        ax31.plot(V[lim:-(lim)], dLIA1[lim:-lim], '--', c='c', label='derivative')
        ax3.set_ylabel('dG/dV', fontsize = 16)
    elif dgdv == 'No':
        ax22.plot(V[:-1],LIA2_int, '--', c='c',label='integrated')
        ax3.plot(V, Iten, c='k')
        ax3.set_ylabel('Intensity [1/V]', fontsize = 16)
        ax31.plot(V[lim:-(lim)], dLIA1[lim:-lim]/LIA1[lim:-(lim)], '--',\
                  c='c', label='derivative')
    ax1.set_ylabel('current [\u00b5A]', fontsize = 16)
    ax1.tick_params(axis='both' , labelsize = 12.0)
    ax2.set_ylabel('dI/dV', fontsize = 16)
    ax2.tick_params(axis='both' , labelsize = 12.0)
    ax22.tick_params(axis='both' , labelsize = 12.0)
    ax3.set_xlabel('Bias [mV]', fontsize = 16)
    ax3.tick_params(axis='both' , labelsize = 12.0)
    ax31.tick_params(axis='both' , labelsize = 12.0)
    if legend == 'yes':
        #ax1.legend()
        ax11.legend()
        ax31.legend()
        ax22.legend()
        #ax2.legend()
    f2.tight_layout()
    #=====================================
    #plt.close('all')
    f3 = plt.figure(33)#, figsize=(6,8.5))
    ax1 = f3.add_subplot(111)
    ax2 = ax1.twinx()
    ax1.plot(V, LIA1, c='k', label='measurement')
    ax2.plot(V[lim:-(lim)], dI[lim:-lim],'--', c='c', label='derivative')
    ax2.set_ylabel('dI/dV [\u00b5S]', fontsize = 16)
    ax1.set_ylabel('G [G0]', fontsize = 16)
    ax1.set_xlabel('Bias [mV]', fontsize = 16)
    if legend == 'yes':
        ax2.legend()
    return f2, f3

def fig_fit(x, y, fi_para, ylab, ti, laloc='lower right',*ax2):
    """
    fi_para : fit parameters
    ylab   : string of the y-label
    ti: string- title
    laloc: label location


    """
    plt.close('all')
    comps, bf, ing = fi_para
    #x = x*1000
    fig_fit = plt.figure(figsize=[6.4, 4.8])
    ax1 = fig_fit.add_subplot(111)
    ax1.plot(x, y,c='k', label='measurement')
    ax1.plot(x, bf, c='r',lw=1, label='best fit')
    ax1.plot(ing[0], ing[1], 'vr')
    #Title = 'Fit'
    #ax1.set_title(Title)
    if ax2:
        ax2 = ax1.twinx()
        for name, comp in comps.items():
            if type(comp) is type(np.array([])):
                ax2.plot(x, comp, '--', lw=0.8, label=name)
    else:
        if 'offset_' in comps.keys():
            ofs = comps['offset_']
            for name, comp in comps.items():
                if type(comp) is type(np.array([])):
                    ax1.plot(x, comp+ofs, '--', lw=0.8, label=name)
                    #ax1.axhline(0, xmin=min(x), xmax=max(x), linestyle=':',lw=0.6)
        else:
            for name, comp in comps.items():
                if type(comp) is type(np.array([])):
                    ax1.plot(x, comp, '--', lw=0.8, label=name)
    lines, labels = ax1.get_legend_handles_labels()
    if ax2:
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc='lower right')
    else:
        ax1.legend(lines, labels, loc= laloc)
    ax1.set_xlabel('Bias [mV]', fontsize = 16)
    ax1.set_ylabel(ylab, fontsize = 16)
    ax1.tick_params(axis='both' , labelsize = 11.0)
    if ti != None:
        plt.title(ti, fontsize = 12)
    fig_fit.tight_layout()
    return fig_fit

def ndr_Plot(x, y, zcr, ylab, mi,facy, xlab='Bias [mV]'):
    """ 
    
        """
    #plt.close('all')
    mi = (mi[0], mi[1]*facy)
    ztxt = float_to_txt(zcr[0])
    fig = plt.figure(25, figsize = [7.68, 5.76] )
    ax11 = fig.add_subplot(111)
    #ax11.plot(x, y, ':' , c='c' , linewidth= 2.2)
    ax11.plot(x, y, c='k' , linewidth= 1.2)
    ax11.scatter(zcr[0], zcr[1], c= 'red', s= 50, marker = "x")
    ax11.scatter(mi[0], mi[1], c= 'darkblue', s= 65, marker = "x")
    ax11.annotate(str(round(mi[1],3)), xy = (mi[0], mi[1]), xytext=(10, -5), \
          color='darkblue', textcoords = 'offset points', ha='left', \
              va='bottom', fontsize = 14)
# =============================================================================
#     ax11.axvline(zcr[0][0], linestyle = '--', c = 'r', linewidth= 0.6)
#     ax11.axvline(zcr[0][1], linestyle = '--', c = 'r', linewidth= 0.6)
# =============================================================================
    ax11.axhline(0, linestyle = '-.', c = 'r', linewidth = 0.6)
    for label, x, y in zip(ztxt, zcr[0],zcr[1]):
        ax11.annotate(label, xy = (x , y), xytext=(10, 5),color='r',\
        textcoords = 'offset points', ha='right', va='bottom',\
        fontsize = 14)
    ax11.set_xlabel(xlab,fontsize=18)
    ax11.set_ylabel(ylab,fontsize=18)
    ax11.tick_params(axis='both' , labelsize = 16.0)
    fig.tight_layout()
    #ax21.legend()
    return fig

def fig1_single(V,y,legend, signal, ylab, xlab='Bias [mV]'):
    """
    signal 'iets' or 'dGdV', I or didv
    """
    V = V
    if type(y) != list:
        LIA1 = y
    else:
        I,LIA1,LIA2, Iten = y[0], y[1], y[2], y[3]
    f1 = plt.figure(1,figsize=[6.4, 4.8])
    ax1 = f1.add_subplot(111)
    #ax4 = ax3.twinx()
    f1.tight_layout()

    if signal == 'I':
        ax1.plot(V,I, c='k')
        ax1.set_ylabel('current [\u00b5A]', fontsize = 16)
    elif signal == 'didv':
        ax1.set_ylabel(ylab, fontsize = 16)
        ax1.plot(V, LIA1, c='k')
    elif signal == 'dgdv':
        ax1.plot(V, LIA2, c='k')#, label='dG/dV [G0]')
        ax1.set_ylabel(ylab, fontsize = 16)
    elif signal == 'iets':
        ax1.plot(V, Iten, c='k', alpha=1)
        ax1.set_ylabel('Intensity [1/V]', fontsize = 16)
    ax1.set_xlabel(xlab, fontsize = 16)
    ax1.tick_params(axis='both' , labelsize = 12.0)
    ax1.legend([legend])
    f1.tight_layout()
    return f1

def fi_single_ti(V, t, y,ylabel,xlabel,fignumber, sec_axis='n', leg=None):
    f1 = plt.figure(fignumber,figsize=[6.4, 4.8])
    ax1 = f1.add_subplot(111)
    f1.tight_layout()
    ax1.plot(t,y, c='k')
    ax1.set_xlabel(xlabel, fontsize = 16)
    ax1.set_ylabel(ylabel, fontsize = 16)
    ax1.tick_params(axis='both' , labelsize = 12.0)
    if sec_axis == 'y':
        ax2 = ax1.twiny()
        ax2.plot(V,y, linestyle='None')
        ax2.set_xlabel('Bias [mV]', fontsize = 16)
        ax2.tick_params(axis='x' , labelsize = 16.0)
    else:
        pass
    ax1.legend([leg])
    f1.tight_layout()
    return f1

def fifi_l(x,y,c_dict,results,sig='Intensity [1/V]', *ax2): #fig fit lin
    """-figure fit linear- ing is initial guess peak positions"""
    plt.close('all')
    #x = x*1000
    fig_fit = plt.figure(10)
    ax1 = fig_fit.add_subplot(111)
    ax1.plot(x, y,c='k', label='measurement')
    fit = c_dict['linear']
    ax1.plot(x, fit, c='r',lw=1, label='fit')
    if ax2:
        ax2 = ax1.twinx()
        for name, comp in c_dict.items():
            if type(comp) is type(np.array([])):
                ax2.plot(x, comp, '--', lw=1, label=name)
    else:
        for name, comp in c_dict.items():
            if type(comp) is type(np.array([])):
                ax1.plot(x, comp, '--', lw=1, label=name)
    lines, labels = ax1.get_legend_handles_labels()
    if ax2:
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax2.legend(lines + lines2, labels + labels2, loc='lower right')
    else:
        ax1.legend(lines, labels, loc='lower right')
    ax1.set_title(['I = ' + str(results[0])+' *Ub' + '+-'+ str(results[1])])
    ax1.set_xlabel('Bias [mV]', fontsize = 16)
    ax1.set_ylabel(sig, fontsize = 16)
    ax1.tick_params(axis='both' , labelsize = 12.0)

    fig_fit.tight_layout()
    return fig_fit

def Coul_Plot(x, y, prep, ylab, xlab='Bias [mV]', txtp =1.0, legpos='best'):
    """
    x is the x array
    y: y array to plot
    prep: dictionary containg the Information
    ylab: wich plot are we looking at?
    xlab: label of x-axis
    txtp: textposition, shifts the string on the x axis left and right.
    ...
    """
    if ylab == 'iets':
        ylab = '$Intensity [V^{-1}]$'
    elif ylab == 'dGdV':
        ylab = '$dG/dV \ [(G/G_{0})*V^{-1}]$'
    elif ylab == 'didv':
        ylab = '$dI/dV \ [G/G_{0}]$'
    elif ylab == 'didmv':
        ylab = '$dI/dV \ [mG/G_{0}]$'

    pd = prep['extrema']
    line = prep['Uline']
    xp = (line[0][-1]-line[1][0])/2
    #yp = min(pd[1])
    pdtxt = float_to_txt(pd[0])
    leg = prep['onoff'][1]
    uxy = prep['U']['U_diagram_xy']
    Ustr = prep['U']['Ustring']
    fig = plt.figure(25, figsize = [7.68, 5.76] )
    ax11 = fig.add_subplot(111)
    ax11.plot(pd[0],pd[1], 'bo')
    #ax11.plot(line[0], line[1] , 'r')
    #ax11.plot(x, y, ':' , c='c' , linewidth= 2.2)
    ax11.plot(x, y, c='k' , linewidth= 1.2)
    ax11.axvline(pd[0][0], linestyle = '--', c = 'r', linewidth= 0.5)
    ax11.axvline(pd[0][1], linestyle = '--', c= 'r', linewidth= 0.5)
    for label, x, y in zip(pdtxt, pd[0],pd[1]):
        ax11.annotate(label, xy = (x , y), xytext=(10, 5), color = 'b', \
        textcoords = 'offset points', ha='right', va='bottom', fontsize = 14)

    ax11.annotate(Ustr, xy = (xp, uxy[1]), xytext = (xp*txtp, uxy[1]*1.01),\
                  color='r', ha='center', va='bottom',fontsize = 14)
    #aro = [min(line[0]),max(line[0])] # arrow
    #for i in range(len(aro)):
    ax11.annotate('', xy = (line[0][0], line[1][0]), xycoords='data', \
          xytext=(line[0][-1], line[1][0]), arrowprops=dict(arrowstyle= '<->',\
                                                    color = 'red',)) #'<|-|>'))
    #ax11.arrow(line[0][0], line[1][0], arlen, 0, color = 'r', head_width = 9E-3, head_length = 2, overhang=-10)
    ax11.set_xlabel(xlab,fontsize=18)
    ax11.set_ylabel(ylab,fontsize=18)
    ax11.tick_params(axis='both' , labelsize = 16.0)
    #fig.tight_layout()
    ax11.set_xlabel(xlab,fontsize=18)
    ax11.set_ylabel(ylab,fontsize=18)
    ax11.tick_params(axis='both' , labelsize = 16.0)
    ax11.legend([leg[0] + '\n' + leg[1]+ '\n' + leg[2]],fontsize = 14,  \
                markerscale=2, loc = legpos)
    fig.tight_layout()
    return fig











