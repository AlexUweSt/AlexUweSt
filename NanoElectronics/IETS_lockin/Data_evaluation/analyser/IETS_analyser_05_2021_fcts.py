# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 10:41:45 2019

@author: strobe46
"""
#import analytic_wfm
import datetime as dt
import os
import scipy.constants as spc
#import scipy.integrate as integ
import matplotlib.pyplot as plt
import numpy as np
import fcts_analyser as fct
import json
#from scipy.fft import fft, ifft
#from scipy.signal import blackman, blackmanharris, hann, tukey


##############################################################################



#%% programm start - paths and datetime =======================================
##########################################################################
Name = 'C60_300ms_0IETS0001_9sdelay_0,2Uspt'

#name = Name
#root = 'C:/Users/Strobel/Nextcloud/Evaluation_data/'
#root = 'C:/Users/Strobe46/ownCloud/Evaluation_data/'
#root = 'C:/Users/Strobe46/ownCloud/Ge_Nanowire_IETS/'
root = '/home/filipk/Documents/Alex/Evaluation_data/IETS_measure/'
path = root + '2019/Oct_Nov_eval_2021/all_measurements/'
#path_add = 'E11/10_11/58/'
#path_add = 'E9/14_10_Coul_Kond/19/'
path_add = ''
ocn = 'OC_0/'
#########################################################################

Date = str(dt.date.today())
name = Name # fct.changename(Name, path_add)
print('new folder name:  ' , name)
# save path

mainpath = path + path_add
ocpath = path + path_add + ocn
#os.rename(ocpath+Name, ocpath+name)
#fpath = ocpath + name
fpath = mainpath + name
#dfpath = fpath + '/processed_arrays_files' # data file path

# Measurement Parameters and constants  =====================================
plp = fct.make_plots_path(fpath, 'analyzer_plots')
G0 = spc.physical_constants['conductance quantum'][0]
curve_attributes, gen_infos, curves, Extrema = {}, {}, {}, {}

#%% load data from measurement arrays ========================================
##########################################################################
cutar = 0 # Cuts Datapoints in the beginning of the array if there are arteffacts (value should be between 0 and 15)
loadfile = 'dict'
##########################################################################

# read from the in python saved files:
if loadfile == 'arr.txt':
    Vall, yall0, Header, sec = fct.load_full_analyze(fpath, name, cutar,\
                                                 1E3, Acamp='auto',ti='no')
elif loadfile == 'dict':
    Vall, yall0, Header, sec = fct.load_processed(fpath, name, cutar, \
                                                  ti='no', onefile = 'true')
# measurement_parameters

Legend,T_insitu,G_OC,Pkdiff  = fct.Parameters(Header, curve_attributes)

#oc_bias = round(np.loadtxt(fpath + '/neghalfsweep.txt')[8][0],2)*Ufac
oc_bias = round(np.loadtxt(fpath + '/processed_arrays_files/1st_sweep.txt')[8][0],2)*1E3

curve_attributes['name'] = name
curves['name'] = name
Extrema['name'] = name
gen_infos['theo_resolution'] =  Pkdiff
gen_infos['eval_Date'] = Date
gen_infos['oc_bias[mV]'] = oc_bias

#%% plot overview diagramm ====================================================
###########################################################################
G_fac = 1E0
###########################################################################
gen_infos['G_scale_fac'] = G_fac # G-scale factor
gen_infos['U_scale_fac'] = 1E3 # U-scale factor
curve_attributes['gen_inf'] = gen_infos

gylab = fct.G_y_label(G_fac)
dgdvlab = '$dG/dV \  [(G/G_{0})V^{-1}]$'

yall = fct.G_scale(yall0, G_fac) # multiply conductance by given factor

f1 = fct.fig1(Vall[0],yall[0],T_insitu,'dgdv', gylab)
#f1.savefig(plp + '/stack_whole_dgdv.png',dpi=200,format='png')

#%% plot ups and downs also ===================================================
###########################################################################
cl = ['full', '1st', '2nd']
udc = ['pos','neg'] # up down curve
###########################################################################
fud = fct.fig_up_down(Vall, yall, ['I','didv','dgdv', 'iets'],udc,cl,gylab)
fud.savefig(plp +'/up_down_sweeps_dgdv.png',dpi=250,format='png')

#%% plot ups and downs stack all ===================================================
###########################################################################
# =============================================================================
# cl = ['full', '1st', '2nd']
# udc = ['pos','neg'] # up down curve
# ###########################################################################
# fud = fct.fig_up_down(Vall, yall, ['I','didv','dgdv', 'iets'],udc,cl,gylab)
# fud.savefig(plp +'/up_down_sweeps_all.png',dpi=250,format='png')
# =============================================================================
#%% Iets scatterplot ======================================================
# iets can also be: [I],[dIdV],[dGdV]
fsc = fct.fig_scatter(Vall, yall, ['iets'], ['pos','neg'])
fsc.savefig(plp+'/scatter_iets_up_down.png',dpi=200,format='png')

#%% scatter dgdv Spectrum ====================================================
fsc_dgdv = fct.fig_scatter(Vall, yall, ['dgdv'], ['pos','neg'])
fsc_dgdv.savefig(plp+'/scatter_dgdv_up_down.png',dpi=200,format='png')

#%% scatter didv
fsc_didv = fct.fig_scatter(Vall, yall, ['didv'], ['pos','neg'])
fsc_didv.savefig(plp+'/scatter_didv_up_down_.png',dpi=200,format='png')
#%% take a closer look.
############################################################################
licl = [-35, 35] # lim closer
############################################################################
wtt = fct.window(Vall[0], yall[0], licl )
zfig = fct.fig1(wtt[0],wtt[1:],T_insitu,'dgdv', gylab)
#zfig.savefig(plp+'/_zoom_in_stack.png',dpi=200,format='png')

#%% attribute entry 2
eye_rating = {}
eye_rating['regime'] =  'iets'
eye_rating['symmetry'] = 'symm'
eye_rating['quality'] = (1,1,2) # like german grades.. i, G, iets
eye_rating['IV_shape'] = 'linear'
eye_rating['G_shape'] = 'Coulomb'
eye_rating['1st-2nd_cons'] = ('no','no') # means up-downsweep consistent for didv and dgdv curve
eye_rating['crack_pos'] = [-155]#[None]
eye_rating['1_2_ef'] = 'yes'
eye_rating['iets_noise'] = 2 # like german grades
eye_rating['noise_thresh'] = (-100, 120)#, (-200,160)
#eye_rating['oscillations'] = 'yes'

#############################################################################
curve_attributes['eye_rating'] = eye_rating # add eye rating to curve attrivutes

#%% window of interests========================================================
#############################################################################
Vlim1 =(-2300, 2300)
############################################################################
pepda = {} # peak detection parameters
misc = {}   
pepda['IETS_window [mV]'] = Vlim1
w = fct.window(Vall[0], yall[0], Vlim1) #[V,V2,I,LIA1,LIA2,Iten]
Vw1, yw  = w[0], w[1:]
Vw = Vw1
#tt = fct.Vtest(Vw) # test if bias is linear..
wn, wp = fct.window(Vall[1],yall[1],Vlim1), fct.window(Vall[2],yall[2],Vlim1)
yw_all, Vw_all = [yw,wn[1:],wp[1:]], [Vw1,wn[0],wp[0]]

#%%= plot stack up_down window =======================================
#fig_check = plt.scatter(Vw,yw[3],s=5)
fud_wi_=fct.fig_up_down(Vw_all, yw_all, ['I','dI/dV','iets'], udc,cl,gylab)
#%% up down dgdv ==============================================
fud_wi=fct.fig_up_down(Vw_all, yw_all, ['I','dI/dV','dgdv'], udc,cl,gylab)
#%% up down stack all
fudstal=fct.fig_up_down(Vw_all,yw_all,['I','didv','dgdv','iets'],udc,cl,gylab) # fig up down stacked all
#%% plot stack fullsw window
fig_st_ie = fct.fig1(Vw1,yw,T_insitu,'iets', gylab)
#%% plot stack mainsweep 
figsal = fct.fig1_all(Vw1,yw,T_insitu, gylab)
if max(Vall[0]) - Vlim1[1] > 30:
    fig_stack = fct.fig1(Vw1,yw,T_insitu,'dgdv', gylab)
    fud_wi.savefig(plp+'/up_down_wind_dgdv.png',dpi=250,format='png')
    fud_wi_.savefig(plp+'/up_down_wind_iets.png',dpi=250,format='png')
    fig_stack.savefig(plp+'/stack_wind_dgdv.png',dpi=250,format='png')
    fig_st_ie.savefig(plp+'/stack_wind_iets.png',dpi=250,format='png')
    figsal.savefig(plp+'/stack_wind_all.png',dpi=250,format='png')

#%% Interpolation
#############################################################################
ipolf = 1 #interpolation factor
smdeg = 'rough' # smotthfactor
#############################################################################
#sarg = np.argsort(Vw) # sort arguments
Vws, yws = fct.sorter(Vw, yw)
Vspl, yspl = fct.interp(Vws,yws,ipolf)
#Vspl, yspl = fct.interp1d(Vw,yw,ipolf) # if other fct is not working "Error on input data"
pepda['interpol_fac'] = ipolf

#= smoothing ================================================================
#fp=[513,6,0,0,'interp'] #window, polyorder, deriv, delta (for deriv)
ysm = fct.smooth(yspl,fct.smooth_param(Vspl,smdeg))
f5 = fct.fig_compare(Vw1, Vspl, yw, ysm, ['as measured','smooth'],smdeg,gylab)
f5.savefig(plp+'/smooth_compare.png',dpi=150,format='png')
pepda['smooth_degree'] = smdeg

#%%  peakdetection process 1 (dgdv) ===================================================
### parameters and energy limiting

#pedet_para_ = [25,25] #[ymax/detet_para[0],len(y)/pedet_para[1]]
#############################################################################
li = [60, 1020] # limit where peaks are detected --> cutout noise etc.
pp = (80,2E-2) # peak detection parameters (lookahead, delta)
#############################################################################
pepda['Peak_det_limit'] = li
pepda['peakf_param_diff_dgdv'] = pp #

# detection and figure
P1 = fct.Peak_diff(Vspl, ysm[2], Pkdiff,pp,fac= 10,lim=li) # output: [peaks,dips], extrema, pos_Extr, [sypecurve,sydicurve]
f2 = fct.peakPlot(Vspl, ysm[2], pd = P1[0], spec ='whole', ylab= dgdvlab)

print('adjust peak positions!')
print('\n adjust peak positions! \n adjust peak positions!')

#%% shifted  spectrum preperation =====================================
###############################################################################
##====  do it again with the shifted potential ==============================
#############################################################################
#Vshift, iets_shift_pd = fct.xy_shift(P1,2,2) #dip,peak count from left
iets_shift_pd = 0 # shift spectrum on y- axis
Vshift = 0.000 # shift spectrum on x-axis
#############################################################################
misc['v_y_shift'] = (Vshift,iets_shift_pd)

if Vshift != 0:
    Vsh = Vall[0]+Vshift
    #ysh = y_shift_0(Vsh,y)[0]
    ysh = fct.y_shifter(yall,iets_shift_pd)

    f3 = fct.fig_compare(Vw1,Vw,yw,yw,['before shift','after shift'],3,gylab)
    f3.savefig(plp+'/compare_spec_corretion.png',dpi=250,format='png')
# interpol 2
    spl2 = fct.interp(np.sort(Vw),yw, ipolf)
    Vspl2, yspl2 = spl2[0], spl2[1]

    #Vspl2, yspl2 = Vw, yw ###  ###     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
# smooth shifted spectra
    ysm2 = fct.smooth(yspl2,fct.smooth_param(Vspl2, smdeg))
# =============================================================================
#     f2_1 = fct.fig1_single(Vspl2,ysm2,T_insitu,'dGdV', 55, \
#                            '$dG/dV [(G/G_{0})/V]$')
#     f2_1.savefig(plp+'/dGdV_spec_corrected.png',dpi=250,format='png')
# =============================================================================

    P1 = fct.Peak_diff(Vspl2, ysm2[2], Pkdiff, pp,lim=li) # output: [peaks,dips], extrema, pos_Extr, [sypecurve,sydicurve]
    f2 = fct.peakPlot(Vspl2, ysm2[2], pd = P1[0], spec ='whole', ylab= dgdvlab)
    f2[0].savefig(plp+'/dgdv_extr_whole.png',dpi=200,format='png')

else:
    ysh = yall[0]
    Vsh = Vall[0]
    Vspl2 = Vspl
    ysm2 = ysm
if 'f2' in locals():
    f2[0].savefig(plp+'/dgdv_extr_whole.png',dpi=200,format='png')

    curves['diff_Extrema_dgdv'] =  [i.tolist() for i in P1[3]]
    Extrema['as_measured_Extrema_dgdv'] = [i.tolist() for i in P1[0]]
    Extrema['peakdifff_Extrema_dgdv'] = P1[2]#.tolist()

w = fct.window(Vsh, ysh, Vlim1) #[V,V2,I,LIA1,LIA2,Iten] ; (-0.27,0.27)
Vw, yw = w[0], w[1:]
curves['shifted_spectra'] = np.vstack((Vsh,ysh)).tolist()

#%% integration and derivation ===============================================
###########################################################################
c_cut = 1
###########################################################################
dyw = fct.derivatives(Vspl2, ysm2)
#dyw = deriv_2(Vspl2, ysm2)
yw_int = fct.Integrals(Vspl2, ysm2)
#f3 = fig_deriv(Vspl, ysm, dyw, yw_int, dgdv='yes', lim = 300)
f4,f5 = fct.fig_deriv(Vspl2, ysm2, dyw, yw_int, dgdv ='No', lim = c_cut)
f4.savefig(plp+'/deriv_integ_stack.png',dpi=100,format='png')
f5.savefig(plp+'/deriv_integ_dIdV.png',dpi=100,format='png')

misc['cutoff_deriv_artefacts'] = [c_cut]
curves['deriv_integ']=([i.tolist() for i in dyw],[i.tolist() for i in yw_int])
#curves['deriv_integ']=(dyw,  yw_int)
#%% Peaks antisymmetry analyse iets ===========================================
###########################################################################
ppd = (40, 1E-1)
limd = [-300,300]
###########################################################################
P_div = fct.Peak_divide(Vspl2,ysm2[3],ppd, lim = limd, sig = 'iets', fac=1) #returns:[peaks1,dips1],[peaks2,dips2],specs

# peak_plot show variable:
f10 = fct.peakPlot(Vspl2, ysm2[3], P_div, 'symmetry', 'iets', 'ansy', 'both')
f10[0].savefig(plp+'/iets_extr_antisymm.png',dpi=250,format='png')
#%% 
f10_1 = fct.peakPlot(Vspl2, ysm2[3], P_div, 'symmetry', 'iets', 'bkg', 'both')
f10_1[0].savefig(plp+'/iets_extr_bkgnd.png',dpi=250,format='png')
#%% Peaks antisymmetry analyse dgdv ============================================
###########################################################################
ppd_ = (40,9E-1)
###########################################################################
P_div_ = fct.Peak_divide(Vspl2,ysm2[2],ppd_, lim=limd, sig = 'dgdv', fac=10)

# peak_plot show variable:
f10_2 = fct.peakPlot(Vspl2,ysm2[2], P_div_, 'symmetry', 'dgdv', 'ansy', 'both')
f10_2[0].savefig(plp+'/dgdv_extr_antisymm.png',dpi=250,format='png')
#%% 
f10_3 = fct.peakPlot(Vspl2,ysm2[2], P_div_, 'symmetry', 'dgdv', 'bkg', 'both')
f10_3[0].savefig(plp+'/dgdv_extr_bkgnd.png',dpi=250,format='png')

print('adjust weather its iets or pcs!')
print('\n adjust weather its iets or pcs! \n adjust weather its iets or pcs!')

#%% peak detection iets whole spec =========================================================
###########################################################################
sign = 'both'
pp_ = (45,9E-1)
###########################################################################
if 'eye_rating' in locals():
    regime = eye_rating['regime']
    
P2 = fct.Peak_diff(Vspl2, ysm2[3], Pkdiff*1000, pp_,fac=10, lim=li)# Pkdiff set to 15 since I wanna see all extrema
# returns: [peaks,dips], extrema, pos_Extr, [sypecurve,sydicurve]

f9 = fct.peakPlot(Vspl2, ysm2[3], P2[0],'whole','iets')#,show=sign)
f9[0].savefig(plp+'/iets_extr_whole.png',dpi=250,format='png')

#%% ad to dictionaries

pepda['peak_param_diff_iets'] = pp_
pepda['lah__del_antisy_iets'] = ppd
pepda['lah__del_antisy_dgdv'] = ppd_

curves.update(P_div_[3])
curves.update(P_div[3])
curves['diff_Extrema_iets'] =  [i.tolist() for i in P2[3]]

Extrema['as_measured_Extrema_iets'] = [i.tolist() for i in P2[0]]
Extrema['peakdiffference_Extrema_iets'] = P2[2]#.tolist()
Extrema['Antisymm_Extrema_dgdv'] = P_div_[1]
Extrema['bkg_subs_Extrema_dgdv'] = P_div_[0]
Extrema['Antisymm_Extrema_iets'] = P_div[1]
Extrema['bkg_subs_Extrema_iets'] = P_div[0]

#%% whole_spectrum with extrema shifted
if Vshift != 0:
    f6_1 = fct.peakPlot(Vspl2 , ysm2[2], P2[0],'whole','dgdv',show=sign)
    f6_1[0].savefig(plp+'/_shifted_spectr.png',dpi=250,format='png')

#%% considered curves and extrema for positive side ===========================
Sycus = fct.symmcurves(P2[3],P2[2],P_div[1], sign) #returns: sig, f, regime, E,Es
#sycu, regime, syExtr, subExtr = Sycus[0], Sycus[2], Sycus[3], Sycus[4]
#fsym = plt.scatter(P2[3][0][0],P2[3][0][1],s=1)
#plt.scatter(sycu[0],sycu[1],s=1.1)
#plt.scatter(syExtr[0],syExtr[1],s=28.1)

#%% window for Au-Dip - peak assuptions =======================================
P3 = P2#Peak_flip(Vspl2, ysm3[3], Pkdiff, peak_params(iets,'middle'))
#P3 returns: [peaks,dips], extrema, pos_Extr, [sypecurve,sydicurve]
yspl3 = yspl2 #y_shift_0(Vspl2,yspl2)[0]
#f12_0 = peakPlot(Vspl2, yspl3[3], P3[0], 'whole', 'iets')
wAu =  fct.window(Vspl2, ysm2, (-22, 26))
VAu = wAu[0]
yAu = wAu[1:]
iets = yAu[3]
dgdv = yAu[2]
f12 = fct.fig1(VAu,yAu,T_insitu,'dgdv',gylab)
#f12.savefig(plp+'/Stack_only_Au_peak.png',dpi=300,format='png')
# fit assumptions ========================================================
#choose what to consider;(all extrema, all dips or all peaks):
xy_pos = P3[1] #P3[0][0], P3[0][1] #--> peaks, dips

# =============================================================================
# #%% fit Gauss =================================================================
# # background options: 'const','linear','quadratic', None
# fit1 = fct.Fit_Gauss(VAu, iets, xy_pos[0], xy_pos[1],'linear')
# # returns: result, comps, best_fit, report, gueses, model, params
# f8 = fct.fig_fit(VAu, iets, fit1[1], fit1[2],6, fit1[4])
# # allpeak fit ==============================================================
# #wfit =  window(Vspl2, ysm3, (-0.021, 0.21))
# #fit1_1 = Fit_Gauss(wfit[0], wfit[4], P3[1][0], P3[1][1],'linear')
# #f8_1 = fig_fit(wfit[0], wfit[4], fit1_1[1], fit1_1[2],8, fit1_1[4])
# #f8_1.savefig(path+'/plots/spec_fit.png',dpi=300,format='png')
# # other fitmodels:
# #%% fit lorentz =============================================================
# fit2 = Fit_lorentz(VAu, iets, xy_pos[0], xy_pos[1],'linear')
# f5 = fig_fit(VAu, iets, fit1[1], fit1[2],5, fit1[4])
# #%% fit split lorentz =======================================================
# fit3 = Fit_slor(VAu, iets, xy_pos[0], xy_pos[1],'linear')
# f7 = fig_fit(VAu, iets, fit3[1], fit3[2],7, fit3[4], 'ax2')
# #%% fit6 voigt ==============================================================
# fit4 = Fit_voigt(VAu, iets, xy_pos[0], xy_pos[1],'linear')# 'const'
# f8 = fig_fit(VAu, iets, fit4[1], fit4[2],8, fit4[4], 'ax2')
# #%% fit skew voigt  =========================================================
# fit5 = Fit_svoigt(VAu, iets, xy_pos[0], xy_pos[1],'linear')
# f9= fig_fit(VAu, iets, fit5[1], fit5[2],9, fit5[4], 'ax2')
# #%%Fit Bose Einstein distribution============================================
# fit6 = Fit_BEinst(VAu, iets, xy_pos[0], xy_pos[1],None)
# f10 = fig_fit(VAu, iets, fit6[1], fit6[2],9, fit6[4], 'ax2')
# #%%Fit DHO==============================================
# fit6 = Fit_Damped_Osc(VAu, iets, xy_pos[0], xy_pos[1],'linear')
# f10 = fig_fit(VAu, iets, fit6[1], fit6[2],9, fit6[4], 'ax2')
# #%%fit doniach ==============================================================
# fit7 = Fit_Doniach(VAu, iets, xy_pos[0], xy_pos[1],None)
# f10 = fig_fit(VAu, iets, fit7[1], fit7[2],9, fit7[4], 'ax2')
# #%%skew gauß fit=============================================================
# fit8 = Fit_sGauss(VAu, iets, xy_pos[0], xy_pos[1],None)
# f13 = fig_fit(VAu, iets, fit8[1], fit8[2],9, fit7[4], 'ax2')
# =============================================================================

#%% Temp calculation if succesfull fit ========================================
#fitp = [0,0]

#Peak_params.append(fitp[1]) # fwhm1+fwhm2/2
#Peak_params.append(fitp[0]) # Temp_fit
# =============================================================================
# fitp = parameters(fit1, ACamp, 'fwhm') # Temp, fwhm, fwhm1, fwhm2
# f9 = fig_fit_T(VAu, iets, fit1[1], fit1[2],6, fit1[4],T_insitu,fitp[0])
# f9.savefig(plp+'/Au_fit_temperature.png',dpi=300,format='png')
# =============================================================================

###############################################################################
###############################################################################

#%% dIdV - analysis

curve_attributes['pepda'] = pepda

print('\n focus now on dIdV spectra')
print('focus now on dIdV spectra')
print('focus now on dIdV spectra')



ped = {}
pecu = {}

if 'Vspl2' not in locals():
    ysh = yall[0]
    Vsh = Vall[0]
    Vspl2 = Vspl
    ysm2 = ysm


#                        differential conductance evaluation

#%%  didv peak detection
###########################################################################
sign = 'both'
ppG = (40,1E-2)
limdi = [-2320,2350] # None #
###########################################################################
PG = fct.Peak_diff(Vspl2, ysm2[1], Pkdiff*1000 , ppG, lim=limdi, fac=100)# Pkdiff set to 15 since I wanna see all extrema
# PG = fct.Peak_diff(Vws, yws[1], Pkdiff*1000 , ppG, lim=limdi, fac=100)
# returns: [peaks,dips], extrema, pos_Extr, [sypecurve,sydicurve]

fpg = fct.peakPlot(Vspl2, ysm2[1], PG[0],'whole',gylab)#,show=sign)
#fpg = fct.peakPlot(Vws, yws[1], PG[0],'whole',gylab)#,show=sign)
fpg[0].savefig(plp+'/didv_extrema.png',dpi=250,format='png')

Extrema['didv_extrema'] = [i.tolist() for i in PG[0]]

#%% Peak Window:
#############################################################################
mpw =  (-140, 40) # multi peak window
pedet_mp = (15, 1E-3) # peak detection multi peak
#############################################################################
m_peaks = {'window' : mpw, 'lah_delta' : pedet_mp}
wmp =  fct.window(Vspl2, ysm2, mpw) # window multi peak
Vmp  = wmp[0] # V zero Bias
ymp = wmp[1:]

Pg_mp = fct.Peak_diff(Vmp, ymp[1],1, pedet_mp, width=[1.,1], lim=None,fac=100) # peak g (conductance) zero Bias
fp_mp = fct.peakPlot(Vmp, ymp[1], Pg_mp[0],'whole','mdidv',show='both')
fp_mp[0].savefig(plp+'/zero_Bias_region.png',dpi=150,format='png')

#%% fit lorentzians =======================================================
wmpf = (-15, 144) #window multi peak fit
###########################################################################
vymp = fct.window(Vspl2, ysm2, wmpf)
xmp,ymp = Pg_mp[0][0][0],Pg_mp[0][0][1]
#xmp,ymp = [6,30,85], [0.081,0.125,-0.006]
#xmp,ymp = Pg_mp[1][0],Pg_mp[1][1]
#fit_mp = fct.Fit_slor(vymp[0], vymp[2], xmp, ymp, 'linear')
fit_mp = fct.Fit_lorentz_si(vymp[0], vymp[2], xmp, ymp, 'const')
#fit_mp = fct.Fit_Gauss(vymp[0], vymp[2], xmp, ymp,'const')# 'linear')
fimp_r = [fit_mp[1], fit_mp[2], fit_mp[4]]

L_csq = fit_mp[0].chisqr # sLo^^rentz_chi_squared
print('\n','chi_square:', L_csq, '\n')
fi3_re = fit_mp[0].fit_report()

fmp = fct.fig_fit(vymp[0], vymp[2], fimp_r, gylab, None, 'best', 'yes') # last yes is fill under curve and is not tested yet 
#%%
fmp.savefig(plp+'/multi_fit.png',dpi=250,format='png')
fmp.savefig(plp+'/multi_fit.svg',format='svg')
#%%
print(fit_mp[0].values.items())
m_peaks['parameters'] = (fit_mp[0].values)

ped['multi_peak_fit'] = m_peaks
pecu['multi_fit'] = fct.dct_np_tolist(fimp_r[0],['best_fit','V-ax'],\
                                [list(fimp_r[1]),list(vymp[0])])


#%% conductances
cnds = {}
Gd0 = fct.Conductance(Vsh,ysh,0.00,G_fac)
Gd50 = fct.Conductance(Vsh,ysh,50,G_fac)
Gd100 = fct.Conductance(Vsh,ysh,100,G_fac)
Gd_m100 = fct.Conductance(Vsh,ysh,-100,G_fac)
Gd_m50 = fct.Conductance(Vsh,ysh,-50,G_fac)
Gd100 = fct.Conductance(Vsh,ysh,100,G_fac)
cnds['Gd0'] = Gd0
cnds['Gd50'] = Gd50
cnds['Gd100'] = Gd100
cnds['Gd_neg_50'] = Gd_m50
cnds['Gd_neg_100'] = Gd_m100

curve_attributes['G_values'] =  cnds


#% Conductance min max and symmerty
# from the didv curve:

cusys = {} # before : curve_symmetries = {}
cusys['dIdV'] = (fct.symm_attr(Vsh, ysh[1],'dIdV_',facx=1E0, \
                                      facy=G_fac))
cusys['I'] = (fct.symm_attr(Vsh, ysh[0], 'I',facx = 1E0, facy=1E0))
cusys['dGdV'] = (fct.symm_attr(Vsh,ysh[2],'dGdV_',facx = 1E0,facy=1E0))
cusys['iets'] = (fct.symm_attr(Vsh, ysh[3],'IETS',facx = 1E0,facy=1E0))


curve_attributes['curve_symmetries'] = cusys

#% Integrals
Intgrls = {}

Intgrls.update(fct.integrator(Vsh,ysh[1],facx=1E3,facy=G_fac, c ='didv'))

Intgrls.update(fct.integrator(Vsh, ysh[2], facx=1E3, facy = 1, c='dgdv'))

Intgrls.update(fct.integrator(Vsh, ysh[3], facx=1E3, facy = 1, c='iets'))

curve_attributes['Integrals'] = Intgrls

print('\n dIdV_general_finished')


#%% NDR analyze


print('\n NDR ANALYZE')

ndrd = {}

#%% defining a Window for NDR peak fitting:
#############################################################################
ndrw =  (-50, 200)
ynd_ = fct.smooth(yspl,fct.smooth_param(Vspl2, 'rough'))
ynd_ = yspl
#############################################################################
wnd =  fct.window(Vspl2, ynd_, ndrw)
Vnd  = wnd[0]
ynd = wnd[1:]
fnd1 = fct.fig1(Vnd,ynd,T_insitu,'dgdv',gylab)
#%% define ndr
###########################################################################
rns = [(50,70), (95,125)]    # ranges
###########################################################################
ndrz = fct.find_y_zero(Vnd,ynd[1],rns)
mini = fct.curve_symmetries['dIdV']['dIdV_min_global']
mini = fct.symm_attr(Vspl2, ysm2[1], 'dIdV_', facx=1E0, \
                     facy=G_fac)['dIdV_min_global']
fndr = fct.ndr_Plot(Vnd,ynd[1], ndrz, gylab,mini,G_fac)
fndr.savefig(plp+'/NDR_analyse.png',dpi=250,format='png')

ndrd['crossings_x_y'] = ndrz
curves['NDR_x_y'], areas = fct.ndr_calc(Vnd, ynd[1], ndrz[0], facy = G_fac)
ndrd['area'], ndrd['area_rang'] = areas

#%%
GIw, GIp = Intgrls['didv_whole_range_'], Intgrls['didv_pos_arm_']
Iwndr, Ipndr = GIw[0] - ndrd['area'], GIp[0] - ndrd['area']

curve_attributes['Integrals']['didv_whole+NDR'] = Iwndr
curve_attributes['Integrals']['didv_parm_+NDR'] = Ipndr

ndrd['GI_rat_whole'] = ((ndrd['area_rang'] *(-1) / Iwndr/(GIw[0]/GIw[1])))
ndrd['GI_rat_parm'] = ((ndrd['area_rang'] *(-1) / Ipndr/(GIp[0]/GIp[1])))

curve_attributes['NDR'] = ndrd
#%% oscilation analyze



print('NDR_finished --> oscilations')



Osc = {}

#%% Oscilations energy and time plots
###########################################################################
cutout = [-199,201]
###########################################################################
[Vcut, secc], yc = fct.IV_cutter([Vsh,sec[:-1]],ysh, cutout)
seccn = secc-min(secc)
# plot
xlab = 'Bias [mV]'
fw15v = fct.fig1_single(Vcut,yc,T_insitu, 'didv', ylab = gylab)
#fcutI = fct.fig1_single(Vcut,yc,T_insitu, 'I', 5)
#%% background ?
print('subtrac "background" signal? \n ####### \n enter "no" or "yes" ')
#%% decide bkg subtract
###########################################################################
subbkg = 'no' #    'no'   # subtract backgroud
fact = 1
###########################################################################
#%% treat curves
if subbkg == 'yes':
    lbkg = fct.subs_lin_bkg(Vcut, yc[1],1E-4,fact)
    yg = yc[1]- lbkg[0] # y minus background
    curve_attributes['Osc_lin_bkg_param'] = lbkg[1]
else:
    yg = yc[1]
yg*1  # y factor to change from mG0 to G0  ????????????????????????????????
fw15v = fct.fig1_single(Vcut,yg, 'didv', ylab = gylab)

print(' \n ---> jump if no oscilation treatment needed')
#%% peak detection GV curve
###########################################################################
ppG = [30, 5E-1]
poscpa = 'middle'
###########################################################################
peakfr = cutout #[10,190] #peak find range
Pg_V = fct.Peak_diff(Vcut, yg, 0.15, ppG,width=[1.8,1], lim = peakfr)
f16 = fct.peakPlot(Vcut,yg,pd=Pg_V[0],spectrum ='whole',ylab='didmv')
#%% delta peak bias dip calc and plot
###########################################################################
ylb='$dG/dV \ [(mG/G_{0})]$' # check if above dgdv lab can't bre used instead
###########################################################################
dvpeaks, dvdips = fct.titofr(Pg_V[0][0][0])[1], titofr(Pg_V[0][1][0])[1]
dvpd = [dvpeaks,dvdips]
aveV = int(round(np.average(np.concatenate((dvpeaks,dvdips))))) #average distance of peaks and dips
ampsv, ave_av = fct.amplitude(Pg_V[0])  #list of amplitude differences, average amplitude difference
fpv_leg = 'average delta: '+ str(aveV)+' mV \n amplitude: '\
    + str(round(ave_av,3)) + ' $mG/G_{0}$'                     # string for legend of figure
f_pV = fct.delta_peak_plot(Vcut,yg,Pg_V[0],dvpd,ylb, xlab,fpv_leg, 'lower left')
f_pV.savefig(plp+'/G_Osc_voltage_.png',dpi=250,format='png')
#%% time osci plot
fw15t = fct.fig1_single(seccn,yg, 'didv',gylab, xlab = 't [sec]')
#%% peakdetect & plot Gt curve
tfr = fct.time_limit(peakfr, Vcut, seccn,yg)[0]
Pg_t = fct.Peak_flip(seccn, yg, 15, ppG,width=[10.8,10])#,lim=tfr)
f16 = fct.peakPlot(seccn, yg, pd = Pg_t[0], spectrum ='whole',\
                   ylab='didmv', xlab = 'time [s]')
#%% calc frequency and plot peaks and dleta peaks of Gt curve
peakf, dtpeaks = fct.titofr(Pg_t[0][0][0])
dipf, dtdips = fct.titofr(Pg_t[0][1][0])
ampst, ave_at = fct.amplitude(Pg_t[0])
dtpd = [dtpeaks,dtdips]
avet = int(round(np.average(np.concatenate((dtpeaks,dtdips)))))
avef = round((1/avet)*1E3)
fpt_l =  'average delta: '+ str(avet)+' s -> '+str(avef) +\
    ' mHz \n amplitude: ' + str(round(ave_at,3)) + ' $G/G_{0}$' #freq plot legend
f_pt = fct.delta_peak_plot(seccn,yg,Pg_t[0],dtpd, ylb,'time [s]',fpt_l,\
                           'lower left')
f_pt.savefig(plp+'/G_Osc_time_.png',dpi=200,format='png')
#%% curve attributes
Osc['Osc_window'] = cutout
Osc['Ocs_pedet_param'] = ppG, poscpa
Osc['Osc_peakfind_range'] = peakfr
Osc['average_osc_peak_dif_v'] = aveV
Osc['average_osc_peak_dif_time'] = avet
Osc['average_osc_peak_amp'] = ave_at # time and voltage should be same
Osc['delta_t_measure'] = fct.titofr(sec)[1]
#curve_attributes['average_osc_peak_amp_bias'] = ave_av

curves['Window_Oscillations_V_G_t']= [Vcut,yg, secc]
#curves['Window_Osc_time']= Pg_t[0]

Extrema['Extrema_Osc_Bias']= Pg_V[0]
Extrema['Extrema_Osc_time']= Pg_t[0]

#%% fourier window (blackman/ hann / etc..  edge treatment)
###########################################################################
ftam = 1E4
fft_fra = [5E-3,0.035]
###########################################################################
tilimarr = fct.time_limit([-200,200],Vall[0],sec[:-1],yarr=yall[0][1])
fft_win,secf,vcf,ycf = tilimarr
f_fft_1 = fct.fi_single_ti([0], secf,ycf, 'Amplitude [a.u.]','time [s]',1,'n')
fft_param = fct.fourier(sec[:-1], yall[0],fft_win, 'hann', fft_fra,ftam)
# output is: refft, refftw, fvec, tc, yc, prod, maxfr
#%% fft plots
###########################################################################
xla, yla = 'f [Hz]','Amplitude [a.u.]'
###########################################################################
refft, refftw, fvec = fft_param[0:3]
sm_fft = fct.smooth(fft_param[0:3], fct.smooth_param(fft_param[0], 'middle'))
#xla = 'freq [\u03bcHz]'
lgd = 'maximum @ ' + str(int((fft_param[6]*1000).round(1))) + ' mHz'
f_fft_2 = fct.fi_single_ti([0],sm_fft[2],sm_fft[0],yla,xla,1,leg=lgd)
f_fft_2.savefig(plp+'/G_Osc_fft_.png',dpi=250,format='png')
#%% fft plots 2 (no window curve, etc less important)
#f_fft_2 = fct.fi_single_ti([0],fvec,refft,'Amplitude',xla,1,'n')
#f_fft_3 = fct.fi_single_ti([0],fvec,refftw,'Amplitude',xla,1,'n')
#%% ad to curves and attributes
fft_pa_nam = ['real_fft','real_fft_win','freq_vec','ttf_time','lia_fft_win']
for i,j in enumerate(fft_param[:-2]):
    curves[fft_pa_nam[i]] = j

#%% dictionary entries
Osc['fft_param_window'] = fft_win
Osc['fft_param_window_function'] = 'hann'
Osc['fft_param_amplitude_amp'] = ftam
Osc['fft_param_freq_range'] = fft_fra
Osc['fft_max_freq_[Hz]'] = fft_param[6]

#%% Osc treatment finisched

curve_attributes['oscillations'] = Osc



print('Osc treatment finisched. \n next: IV-curve fit -roughly ')




#%% fitting the IV_curve (for now only just a few curves possible)
IV_treatm = {} #-->   IV - treatment
ivfit = fct.Fit_linear(Vall[0], yall[0][0])
ivfr=round(ivfit[6]['slope'].value*1E0,9),round(ivfit[6]['intercept'].value,9) # IV fit results
print(ivfr)
f_xy = fct.fifi_l(Vall[0], yall[0][0],ivfit[1],ivfr, 'current [\u00b5A]')

IV_treatm['IV_fit_param_[1V_A]'] = ivfr[0]/1E3,ivfr[0]/1E6

#%% IV curve fit

curve_attributes['IV_treatment'] = IV_treatm


print('IV-curve fit finisched. \n next: Coulomb blockade analyze ')





#%%


print('Coulomb_blockade_analyze')


cold = {}  # coulomb blockade dictionary

#%% Coulom blockade evaluation
###########################################################################
cutout = [-100,160]
xlab = 'Bias [mV]'
###########################################################################
wco = fct.window(Vspl2, ysm, cutout)
Vcut = wco[0]
yc = wco[1:]
# plot

fw15v = fct.fig1_single(Vcut,yc,T_insitu, 'didv', ylab = gylab)
#%% plot coul diagramm
###########################################################################
ppCo = (5,1E-1)
###########################################################################
G_dict = fct.symm_attr(Vcut, yc[1], 'dIdV_') # check if something here can be turnded into a variable

Pg_Co = fct.Peak_diff(Vcut, yc[1], 1, ppCo, width=[1.,1], lim= None)

Pco = fct.peakPlot(Vcut, yc[1], pd = Pg_Co[0], spec ='whole', ylab= 'didv')
#%% create coulomb parameters
###########################################################################
pick_ex = [[0,1],[0]] # pick extrema for coulomb eval [[dips], [peaks]], count from left
###########################################################################
#fCo = fct.peakPlot(Vcut, yc[1],pd=Pg_Co[0],spectrum ='whole',ylab='didmv')

pcoul = fct.prep_Coul_(Pg_Co[0], G_dict, pick_ex)

fcoul = fct.Coul_Plot(Vcut,yc[1],pcoul,ylab='didv',xlab='Bias [mV]',txtp=0.5,\
                     legpos = 'center')


#%% curve_attributes coulomb
fcoul.savefig(plp+'/Coulomb_analyze.png',dpi=300,format='png')
cold['U_coulomb'] = (pcoul['U']['U_coulomb'])
cold['coulomb_peaks&valley'] = (pcoul['extrema'])
cold['Coulomb_on_off'] = (pcoul['onoff'])
cold['extrema_ind'] = pick_ex
cold['lah_delta_coul'] = ppCo
curve_attributes['coulomb'] = cold

#%% Kondo analyze



print ('\n coulomb_analyse finisched! \n')
print('Coulomb analyze finisched --> Kondo analyze   \n Kondo \n Kondo')



kondd = {}2
kocu = {}

#%% defining a Window for kondo peak fitting:
#############################################################################
konw =  (-65, 65)
#yKo_ = fct.smooth(yspl,fct.smooth_param(Vspl2, 'fine'))
#############################################################################
wKo =  fct.window(Vspl2, ysm2, konw)
VKo  = wKo[0]
yKo = wKo[1:]
f12 = fct.fig1(VKo,yKo,T_insitu,'dgdv',gylab)

#%% Kondo Fitting

print('choose fct to fit the Kondo Peak (Lorentz is theoretical  best!')


#%% fit split lorentz =======================================================
###########################################################################
wkl = (-20, 20) #window kondo (split) lorentz
###########################################################################
wKol = fct.window(Vspl2, ysm2, wkl)
ped['window_Kondo_lorenz'] = wkl
zbfit = {'window': wkl}
fit3 = fct.Fit_lorentz_si(wKol[0], wKol[2],  xmp, ymp, 'linear')
#fit3 = fct.Fit_slor(wKol[0], wKol[2],  xmp, ymp, 'linear')
fi3_r = [fit3[1], fit3[2], fit3[4]]
L_csq = fit3[0].chisqr # sLo^^rentz_chi_squared
print('\n','chi_square:', L_csq, '\n')
fi3_re = fit3[0].fit_report()
f7 = fct.fig_fit(wKol[0], wKol[2], fi3_r ,gylab,None,'best')
print(fit3[0].values.items())
#%% fit Fano normalized ==============================================================
###########################################################################
wkfi =   (-35, 22)
###########################################################################
xy_pos = [np.array([0]), np.array([Gd0])]
wKofi = fct.window(Vspl2, ysm2, wkfi)  #
#wKofi = fct.window(Vspl2, yKo_, wkfi)
kondd['window_Fano_norm'] = wkfi
fitFn = fct.Fit_Fano_n(wKofi[0],wKofi[2],xy_pos[0],xy_pos[1],'linear',20)# 'const'
cen = fitFn[0].values['fan_1center']
ceyn = fitFn[2][np.where(np.isclose(wKofi[0],cen, atol=5E-1))[0][0]]
fiFn_r = [fitFn[1], fitFn[2], [cen, ceyn]]
Fan_csq = fitFn[0].chisqr # Fano _normalized_chi_squared
print('\n','chi_square:', Fan_csq, '\n')
f8 = fct.fig_fit(wKofi[0], wKofi[2], fiFn_r, gylab,None, 'best')
print(fitFn[0].values.items())
#%% fit Fano ==============================================================
###########################################################################
#wkf = (-42, 17)
###########################################################################
wKof = wKofi # fct.window(Vspl2, yKo_, wkf)
#kondd['Kondo_window'] = wkfi
fitF = fct.Fit_Fano(wKof[0], wKof[2], xy_pos[0], xy_pos[1],'linear',50)# 'const'
ce = fitF[0].values['fa_1center']
cey = fitF[2][np.where(np.isclose(wKof[0],ce, atol=6E-1))[0][0]]
fiF_r = [fitF[1], fitF[2], [ce, cey]]
Fa_csq = fitF[0].chisqr # Fano_chi_squared
print('\n','chi_square:', Fa_csq, '\n')
f8 = fct.fig_fit(wKof[0], wKof[2], fiF_r, gylab, None, 'best')
print(fitF[0].values.items())


#%% evaluate fit split  Lorentz
peti = fct.Peak_fig_title(fit3[0].values, 'lz') #peak figure title

Ko_ff = fct.fig_fit(wKol[0], wKol[2], fi3_r ,gylab,peti,'best')# Kondo final figure
#Ko_ff.savefig(plp+'/fit_peak_LORENTZ.png',dpi=300,format='png')
zbfit['parameters'] = (fit3[0].values)
ped['peak_0V'] = zbfit
pecu['peak_0V'] = fct.dct_np_tolist(fi3_r[0],['best_fit','V-ax'],\
                              [list(fi3_r[1]),list(wKol[0])])

#%% evaluate fit - Fano_normalized
fitKo_r, Ko_fwhm = fitFn[0].best_fit, fitFn[0].values['fan_1sigma']
KoTemp_ = round((spc.e * Ko_fwhm/1000 )/(2*spc.k * T_insitu),1)

KoTemp = round((spc.e * Ko_fwhm/1000 )/(2*spc.k * T_insitu),1)
Koti = fct.Kondo_fig_title(fitFn[0].values, KoTemp, T_insitu, 'fano_n-')

Ko_ff = fct.fig_fit(wKofi[0], wKofi[2], fiFn_r , gylab, Koti, 'best')# Kondo final figure
Ko_ff.savefig(plp+'/Kondo_peak-FANO_norm.png',dpi=250,format='png')

kondd['parameters_Fano_norm'] = (fitFn[0].values)
kondd['Kondo_temp_Fano_norm'] = (KoTemp)

kocu['Fano_normalized'] = fct.dct_np_tolist(fiFn_r[0],['best_fit'],\
                                      [list(fiFn_r[1])])

#%% print only fitte fano_n curve
fx8 = fct.fig1_single(wKofi[0],fiFn_r[0]['fan_1'],'fano_fit_n','didv',gylab)
fx8.savefig(plp+'/Kondo_peak_FANO_norm_fit.png',dpi=120,format='png')

#%% evaluate fit Fano
fitKo_r, Fa_fwhm = fitF[0].best_fit, fitF[0].values['fa_1sigma']
KoTemp = round((spc.e * Fa_fwhm/1000 )/(2*spc.k * T_insitu),1)
Koti = fct.Kondo_fig_title(fitF[0].values, KoTemp, T_insitu, 'fano')

Ko_ff = fct.fig_fit(wKof[0], wKof[2], fiF_r ,gylab, Koti,'best')# Kondo final figure
Ko_ff.savefig(plp+'/Kondo_peak_FANO.png',dpi=250,format='png')

kondd['parameters_Fano'] = (fitF[0].values)
kondd['Kondo_temp_Fano'] = (KoTemp)

kocu['Fano'] = fct.dct_np_tolist(fiF_r[0],['best_fit'],[list(fiF_r[1])])
#%% print only fano fit
fx8 = fct.fig1_single(wKofi[0],fiF_r[0]['fa_1'],'fano_fit','didv',gylab)
#fx8.savefig(plp+'/Kondo_peak_FANO_bestfit.png',dpi=120,format='png')
#%% end of kondo analyse

curve_attributes['Kondo'] = kondd
curves['Kondo'] = kocu
curves['dIdV_peaks'] = pecu

print('End of Kondo analyse')




#%% finalize


print ('finalyze curve evaluation ')

curve_attributes['misc'] = misc # miscellaneous
curve_attributes['dIdV_peaks'] = ped

####### ##########



#%% save extrema and shifted curves ===========================================

spaths = [mainpath, ocpath, path, fpath]
dictpaths = ['extrema_analysed', 'curves_analysed', 'curve_attributes']
mpa = fct.makepath_local(spaths, name, dictpaths)


#%% shift params to file ======================================================
dicts = [Extrema, curves, curve_attributes]
fct.save_dictionaries(dicts, spaths, name, dictpaths)



#%% END END END END END END


print('__END__ \n   __END__')

#%%
jdna = fpath + '/curve_attributes/' + 'E9_07_10_9___8Hz_sample_freqIETS0012dia.json'
with open(jdna, 'r') as myfile:
    data=myfile.read()
cattr = json.loads(data)#,allow_pickle=True)
