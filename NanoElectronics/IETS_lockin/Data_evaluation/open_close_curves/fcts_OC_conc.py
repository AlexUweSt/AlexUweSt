# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 17:38:36 2021

@author: strobe46
"""

import numpy as np
import os
import matplotlib.pyplot as plt
import fcts_OC_conc_ as fct
import matplotlib as mpl


def cols():
    return ['r','c','g' ,'b','y','lime', 'cadetblue','m','lightgreen', \
      'dimgrey', 'goldenrod', 'slateblue', 'gold', 'deeppink', \
          'steelblue','k', 'orangered', 'magenta', 'maroon', 'saddlebrown']


def conc_OC(files, sw_files):
    tcnt, ic = 0, 0
    tcheck, Gbr, ibr, ic_len = [], [], [], [0]
    swt, swG, ineg, NDRs, cun = extract_swepps(sw_files)
    for count , dat  in enumerate(files):
        #if count in used_OC and os.path.getsize(j[count])>1500:
        OC_file = np.loadtxt(dat , delimiter = '\t', dtype=np.str)
        OC_file = np.char.replace(OC_file, ',' , '.')
        OC_file = np.char.replace(OC_file, 'b', '').astype(np.float64)
        oc = np.transpose(OC_file)
        ts = oc[0]
        #MposOC = oc[1]/1E+7
        vel = oc[2] # velocity of the Motor
        Gval = oc[3]
        #Rval = oc[4]
        LIAG = oc[7]
        tempe = oc[11]
        #LIA1Dmm = oc[5]
        #LIA2Dmm = oc[6]
        plt.close('all')
        #path_ = directory + 'OC_' + str(count) + '/'
        tcnt = tcnt + max(ts)
        ic = ic + len(Gval)
        ic_len.append(ic)
        tcheck.append(tcnt)
        #tiarr.append(ts)
        if count == 0:           #maybe switch if / else again
            ibr.append(ic)
            Gbr.append(Gval[-1])
            tconc, Gvalc, LIAGc, tempec, velc = ts, Gval, LIAG, tempe, vel
        else:
            ibr.append(ic)
            Gbr.append(Gval[-1]) #Gbr --> G breaks --> where an open/close ends
            ts = ts + tcnt
            tconc, Gvalc  = np.hstack([tconc,ts]), np.hstack([Gvalc, Gval])
            LIAGc, tempec = np.hstack([LIAGc,LIAG]), np.hstack([tempec,tempe])
            velc = np.hstack([velc,vel])
            #tempec = np.hstack([tempec,tempe])
    tvec = np.linspace(0, max(tconc)/1000,len(tconc))
    t_br = []
    for i, n in enumerate(ibr):
        t_br.append(tvec[n-1])
    swti = []
    for i, m in enumerate(swt):
        swti.append(tvec[m])
    longarr = [tvec, Gvalc, LIAGc, velc, tempec]
    vals = [longarr,[swti,swG],[t_br,Gbr], ineg,NDRs,cun]
    dct = {'time_vector': tvec, 'Gdrop_conc': Gvalc, 'G_lockin': LIAGc, \
           'velocity':velc,'sweeptime': swti,'sweep_G': swG,'t_oc_end': t_br,\
               'G_oc_end': Gbr, 'datapoints_each_oc': ic_len, \
                   'pos_neg_sweeps': ineg, 'NDR_number_oc': NDRs, \
                       'curve_numbers': cun, 'Temperature':tempec}
    return dct, vals #tvec, Gvalc, LIAGc, velc, swti, swG, t_br, Gbr, tcheck

def extract_swepps(files):
    swt, swG, ineg, NDRs, cun = [], [], [], [], [] # cun -> curve number
    for idx , dat  in enumerate(files):
        if os.path.getsize(dat) < 1.2E3:
            pass
        else:
            V, y, Header, sec = fct.load_main_full(dat, 0, \
                                                   Acamp = 'auto', ti='no')
            param = np.concatenate(([Header[1].astype(np.float)] ,\
                                    [Header[3].astype(np.float)]), axis = 0)
            sweepdat = [param[0,1], int(param[1,0])]
            swG.append(sweepdat[0])
            swt.append(sweepdat[1])

            ####### check for NDR's
            if min(y[1]) > 0:
                sign = 'p'
            elif min(y[1]) < 0:
                sign = 'n'
                NDRs.append((dat[-3:],int(Header[3][1])))
            else:
                sign = 'z'
                #pcurv.append((dat[-3:],int(Header[3][1])))
            ineg.append(sign)
            cun.append(dat[-3:])
            #ipn = np.array(ineg)
    return swt, swG, ineg, NDRs, cun

def choose_OCs(dct, Ocs):
    """
    OCs must to come after each other --> OCs = [1,3], then OC 1, OC 2 and
    OC 3 will be considered, first OC is OC 0.
    """
    sta = dct['datapoints_each_oc'][Ocs[0]]  #start data points
    fi = dct['datapoints_each_oc'][Ocs[1]+1]    # final data points
    t = np.array(dct['time_vector'])
    r = np.arange(0,len(t),1)
    m = (r >= sta) & (r <= fi)
    G = np.array(dct['Gdrop_conc'])[m]
    Gl = np.array(dct['G_lockin'])[m]
    v = np.array(dct['velocity'])[m]
    t = np.array(dct['time_vector'])[m]
    T = np.array(dct['Temperature'])[m]

    sti = np.array(dct['sweeptime'])
    mII = (sti >= min(t)) & (sti <= max(t))
    sti = sti[mII]
    swG = np.array(dct['sweep_G'])[mII]
    swpn = np.array(dct['pos_neg_sweeps'])[mII]
    cuns = np.array(dct['curve_numbers'])[mII]
    cuns = [int(i) for i in cuns]
    cuns = np.array([str(i) for i in cuns])

    tb = np.array(dct['t_oc_end'])
    mIII = (tb >= min(t)) & (tb <= max(t))
    tb = tb[mIII]
    gb = np.array(dct['G_oc_end'])[mIII]
    dct = {'time_vector': t, 'Gdrop_conc': G, 'G_lockin': Gl, 'velocity':v, \
       'sweeptime': sti, 'sweep_G': swG, 't_oc_end': tb, 'G_oc_end': gb, \
           'pos_neg_sweeps': swpn, 'curve_numbers': cuns, 'Temperature': T}
    vals = [[t,G,Gl,v,T],[sti,swG],[tb,gb], swpn, cuns]
    return dct, vals

def choose_time(dct, tim):
    """
    OCs must to come after each other --> OCs = [1,3], then OC 1, OC 2 and
    OC 3 will be considered, first OC is OC 0.
    """
    sta = tim[0]  #start data points
    fi =  tim[1]   # final data points
    t = np.array(dct['time_vector'])
    r = np.arange(0,len(t),1)
    m = (t >= sta) & (t <= fi)
    t = t[m]
    G = np.array(dct['Gdrop_conc'])[m]
    Gl = np.array(dct['G_lockin'])[m]
    v = np.array(dct['velocity'])[m]
    T = np.array(dct['Temperature'])[m]

    sti = np.array(dct['sweeptime'])
    mII = (sti >= min(t)) & (sti <= max(t))
    sti = sti[mII]
    swG = np.array(dct['sweep_G'])[mII]
    swpn = np.array(dct['pos_neg_sweeps'])[mII]
    cuns = np.array(dct['curve_numbers'])[mII]
    cuns = [int(i) for i in cuns]
    cuns = np.array([str(i) for i in cuns])

    tb = np.array(dct['t_oc_end'])
    mIII = (tb >= min(t)) & (tb <= max(t))
    tb = tb[mIII]
    gb = np.array(dct['G_oc_end'])[mIII]
    dct = {'time_vector': t, 'Gdrop_conc': G, 'G_lockin': Gl, 'velocity':v, \
       'sweeptime': sti, 'sweep_G': swG, 't_oc_end': tb, 'G_oc_end': gb, \
           'pos_neg_sweeps': swpn, 'curve_numbers': cuns}
    vals = [[t,G,Gl,v, T],[sti,swG],[tb,gb], swpn, cuns]
    return dct, vals








### figures #####################################################################
def f_conc_oc(concarr,tGswe,tGbr,psw,loca,swpn,cun,si=32,br='axv',fs=(12,5.8),\
              tp = (0,-1.8)):
    """
    concarr: 5 lists with concatenated t, G, Glia, velocity, temp

    tGswe: 2 Lists with the time and G position of the taken sweeps

    tGbr: t and G positions where the OC cycles beginn and end

    psw stands for plot sweep if yes the sweep point are plotted.

    swpn is the array to identify sweep with just positive or also negative
    (NDR) values.

    loca: location of legend

    cun: curve numbers

    size: size of the sweep marker dots

    br: if OC cycles are seperated with x-scatters ('scat')
        or vertical lines ('axv')
    tb: tuple --> position of the sweep nr txt

    """

    t,G,Glia,vel, temp = concarr
    #loca = loca
    #Gval, LIAG = [], []
    tbre, Gbre = tGbr
    tswe, Gswe, cun = np.array(tGswe[0]), np.array(tGswe[1]), np.array(cun)
    mGI = (G > 0)
    mGII = (G < 0)
    mLI = (Glia > 0)
    mLII = (Glia < 0)

    tvp = t[mGI] # time vector positive
    Gvp = G[mGI]

    tvn = t[mGII] # time vector negative
    Gvn = G[mGII]*(-1)

    tlp = t[mLI] # time lia positive
    Glp = Glia[mLI]

    tln = t[mLII] # time lia negative
    Gln = Glia[mLII]*(-1)

    plt.close('all')
   #################################################
    fig, ax = plt.subplots(figsize=fs)

    if min(Gvp) > min(Glp):
        plt.ylim(bottom = min(Glp)- min(Glp)/2)
    else:
        plt.ylim(bottom = min(Gvp) -min(Gvp)/2)
    if max(Gvp) > max(Glp):
        plt.ylim(top = max(Gvp)+max(Gvp)/2)
    else:
        plt.ylim(top = max(Glp)+max(Glp)/2)

    Gp = ax.scatter(tvp, Gvp, s=0.11, color='b',label = 'G')
    Gpd = ax.scatter(tlp,Glp,s=0.11, color='m',label='diff G')

    ax.scatter(tvn, Gvn, s=15.15, color='c')
    ax.scatter(tln,Gln,s=15.15, color='lime')
    if br == 'axv':
        for ii, nn in enumerate(tbre):
            ax.axvline(nn, linestyle = '--', c = 'g', linewidth= 0.6)
    elif br == 'scat':
        ax.scatter(tbre,Gbre,s=200, color='g', marker = "x") # marker = "|"
    if psw == 'yes':
        mp = [i == 'p' for i in swpn]
        mn = [i == 'n' for i in swpn]
        if type(cun) == np.ndarray:
            for label, x, y in zip(cun[mp], tswe[mp], Gswe[mp]):
                ax.annotate(label, xy = (x , y), xytext=(tp),color='k',\
                textcoords = 'offset points', ha='center', va='bottom',\
                fontsize = round(si/10,1))
            for label, x, y in zip(cun[mn], tswe[mn], Gswe[mn]):
                ax.annotate(label, xy = (x , y), xytext = (tp), color='w',\
                textcoords = 'offset points', ha='center', va='bottom',\
                fontsize = round(si/10,1))
        ax.scatter(tswe[mp],Gswe[mp],s=si, color='r', marker = "o")#, label = 'G')
        ax.scatter(tswe[mn],Gswe[mn],si, 'k', "o", label ='NDR')
        swn = ax.scatter(tswe[mn],Gswe[mn],0.15, 'k', "o", label ='NDR')
        plt.legend([Gp, Gpd, swn], ['static G', 'differential G','NDR'], \
           fontsize = 16,loc= loca)#, prop={'si': 6})
    else:
            plt.legend([Gp, Gpd], ['static G', 'differential G'], \
               fontsize = 16,loc= loca)#, prop={'si': 6})
    #plt.ylim(bottom = 1E-9)
    ax.set_yscale('log')
    plt.grid(True, which='major',axis='y',color='g',linewidth=0.2)
    ax.set_xlabel('time [1E3 s]' , fontsize = 20)
    ax.set_ylabel('conductance $[G/G_0]$' , fontsize = 20)
    ax.tick_params(axis='both' , labelsize = 18.0)
    ax2 = ax.twinx()
    ax2.scatter(t , vel , s=0.05 , color = 'k' , alpha=0.01)
    params = {'legend.fontsize': 16, 'legend.handlelength': 1, \
              'legend.markerscale':15}
    plt.rcParams.update(params)
    ax2.tick_params(axis='y' , labelsize = 16.0)
    ax2.set_ylabel('Motorspeed [a.u.]' , fontsize = 16)
    # =============================================================================
    # if abs(MposOC[0]) > abs(MposOC[len(MposOC)-1]):
    #             ax.invert_xaxis()
    # =============================================================================
    plt.tight_layout()
    #plt.show()
    return fig

#def G_temperature (x,y,param, rang, title, design='bone'):
def G_temperature (arr,br,ti=None,dsgn='plasma',lo='best',tlim=None,lab='Gd'):
    """
    arr: 5 lists with concatenated t, G, Glia, velocity, temp
    ti: title string
    lab = select if differetial G ('Gd') or static G should be plottet ('Gs')

    ..
    """
    x, param = arr[0], arr[4]
    if lab == 'Gd':
        y = arr[2]
        legend = 'differetial Conductance'
        yla = '$G \/\ [G/G_0]$'
        Gl = [abs(i) for i in y]
    elif lab == 'Gs':
        y = arr[1]
        legend = 'static conductance'
        yla = '$G \/\ [G/G_0]$'
        Gl = [abs(i) for i in y]
    else:
        y = list(np.array(arr[2])-np.array(arr[1]))
        legend = 'difference'
        yla= lab
        Gl = y# [abs(i) for i in y]


    plt.close('all')
    if tlim == None:
        norm = mpl.colors.Normalize(vmin = min(param), vmax = max(param)+2)
    else:
        norm = mpl.colors.Normalize(vmin = tlim[0], vmax = tlim[1])

    f1 = plt.figure(3,figsize=(11.6, 5.8))
    ax1 = f1.add_subplot(111)
    plt.ylim(bottom = min(y)-min(y)/2)
    plt.ylim(top=max(y)+ max(y)/2)
    if lab == 'Gd' or lab == 'Gs':
        ax1.set_yscale('log')
    f1.tight_layout()
    for ii, nn in enumerate(br):
        ax1.axvline(nn, linestyle = '--', c = 'g', linewidth= 0.4)
    plot = ax1.scatter(x,Gl,c=param,cmap=dsgn,norm=norm ,s=0.6, alpha=0.95)
    #ax1.pcolor([V,param])
    cbar = f1.colorbar(plot)
    cbar.set_label('Temperature [K]', size = 18)
    plt.legend([plot], [legend],  fontsize = 12, loc= lo, markerscale=10)#, prop={'size': 6})
    f1.tight_layout()
    ax1.set_title(ti,size=16)
    #f1.subtitle('OV Bias Conductivity < 0.05 G0',fontsize=16)
    ax1.set_ylabel(yla, fontsize = 20)
    ax1.set_xlabel('time [ks]', fontsize = 20)
    ax1.tick_params(axis='both' , labelsize = 18.0)
    #ax1.legend(['N = ' + str(len(y))+ ' points'])
    f1.tight_layout()
    return f1

def G_special (arr,br,ti=None,dsgn='plasma',lo='best',tlim=None,lab='Gd'):
    """
    arr: 5 lists with concatenated t, G, Glia, velocity, temp
    ti: title string
    lab = select if differetial G ('Gd') or static G should be plottet ('Gs')

    ..
    """
    x, param = arr[0], arr[4]
    if lab == 'Gd':
        y = arr[2]
        legend = 'differetial Conductance'
    elif lab == 'Gs':
        y = arr[1]
        legend = 'static conductance'
    else:
        y = arr[2]-arr[1]
        legend = lab

    Gl = [abs(i) for i in y]

    plt.close('all')
    if tlim == None:
        norm = mpl.colors.Normalize(vmin = min(param), vmax = max(param)+2)
    else:
        norm = mpl.colors.Normalize(vmin = tlim[0], vmax = tlim[1])

    f1 = plt.figure(3,figsize=(11.6, 5.8))
    ax1 = f1.add_subplot(111)
    plt.ylim(bottom = min(y)-min(y)/2)
    plt.ylim(top=max(y)+ max(y)/2)
    ax1.set_yscale('log')
    f1.tight_layout()
    for ii, nn in enumerate(br):
        ax1.axvline(nn, linestyle = '--', c = 'g', linewidth= 0.4)
    plot = ax1.scatter(x,Gl,c=param,cmap=dsgn,norm=norm ,s=0.6, alpha=0.95)
    #plot = ax1.scatter(x,Gl,s=0.6, alpha=0.95)
    #ax1.pcolor([V,param])
    cbar = f1.colorbar(plot)
    cbar.set_label('Temperature [K]', size = 18)
    plt.legend([plot], [legend],  fontsize = 12, loc= lo, markerscale=10)#, prop={'size': 6})
    f1.tight_layout()
    ax1.set_title(ti,size=16)
    #f1.subtitle('OV Bias Conductivity < 0.05 G0',fontsize=16)
    ax1.set_ylabel('$G \/\ [G/G_0]$', fontsize = 20)
    ax1.set_xlabel('time [ks]', fontsize = 20)
    ax1.tick_params(axis='both' , labelsize = 18.0)
    #ax1.legend(['N = ' + str(len(y))+ ' points'])
    f1.tight_layout()
    return f1

