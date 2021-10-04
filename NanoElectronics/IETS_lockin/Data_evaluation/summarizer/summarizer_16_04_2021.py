# -*- coding: utf-8 -*-
"""
Created on Thu Dec 10 10:22:19 2020

@author: strobe46
"""


import numpy as np
import fcts_summarizer as fcts
import datetime

Date = str(datetime.date.today())



#############################################################################
#%%
#== start: directory and load all files within the curves restriction  ======

##############################################################################
# =============================================================================
# root = 'C:/Users/Strobel/Nextcloud/Evaluation_data/'+\
#             'IETS_measurements_revised/2020_measurements/S3_2020/C60/'
# =============================================================================
root = '/home/filipk/Documents/Alex/Evaluation_data/IETS_measure/' + \
           '2020/S3_2020/C60/'
drc = root + '/E11/11_11/59'
directory =  drc #+ '/OC_6'
fcts.fp_exists(drc)
#%% Filter by Curve Number 
curs = [84, 85, 86, 87, 88, 89, 95, 99, 100, 103 ] # if selected by curves
curs = np.arange(51,101,1).tolist()
#del curs[85:88]
#curs.remove(99)
#curs = dinam
##############################################################################
curna = fcts.curname(curs) #string of min max curves
#curna = 'give it another name..'
# load evaluated curves from json files
ldcts = fcts.filt_curv_eval(directory, curs) # list dictionaries --> all files in curs list are loaded 

#%%  # filter by main feature of the experiment ('Kondo), 'Coulomb', 'IETS', 
##############################################################################
mk = ['tests'+curna]  # main feature Key for titles etc
mk_ = ['G_values']
fdct_ = 0 # filter dictionary (the dict to look for mk_)
##############################################################################
# fbk retuns: [atributes, curves, extrema], names    
fids_, nams_ = fcts.fbk(ldcts, key=mk_, fdi=fdct_, inex = None)  # fbk --> filter by key

#%% filter by Conductance  (what </> value) 
##############################################################################
what = 'G0' # choose : GI, G0, GIn, GIp, G50, G100, Gm50, Gm100, Gma
comp = '>' # comparator
val = 0.1 # value 
##############################################################################
if val != None:
    fids, nams = fcts.fbv(fids_, what, comp, val)
else:
    fids, nams = fids_, nams_
#nad['swnr'] = (nams, 'sweep_nr', 'sweep_number')
nad = fcts.new_dict1(fids, nams, mk_)  # new attribute dictionary
#%% choose the Spectra that should be plotted and one parameter (defined by atk
##############################################################################
spec = 'didv' #'peak_0Vbf' #'Kondo_meas'  # spectrum to extract from the filterted curve dictionary 
##############################################################################
evar = fcts.curve_arrays(fids[1], fids[0], spec) #evaluation array, extended x and y lists
       
#arrxx = cu_extract(fids[1], 'shifted_spectra')

#%% cut curves for 2D plots
##############################################################################
cut = [-800,800] #[-160,250] # limits the voltage array of curves
cst = str(min(cut)) + '-' + str(max(cut)) # cut string
nm = 'ma'  # normalisation method
##############################################################################
verts, hiG, loG, Vlim, cutar = fcts.prepare_curves(evar, cut, None)
vertsn, hiGn, loGn, Vlim, cutarn =fcts.prepare_curves(evar, cut,nm) #!!! VisibleDeprecationWarning: Creating an ndarray from ragged nested sequences (which is a list-or-tuple of lists-or-tuples-or ndarrays with different lengths or shapes) is deprecated. If you meant to do this, you must specify 'dtype=object' when creating the ndarray
# ============== normalisation methods: 'ma' --> ; mag -->           ==========
#=============================================================================

    
#%% 


print('plot spectra - multicurve - plot')



#=============================================================================

#%% figure multi plot _ normalized  whole spec
##############################################################################
Tit = mk[0]
xla = 'Bias [mV]'
yla = '$G_{norm}$' #'dG/dV [a.u.]' #  
hrzl = []#[-10,195]#[-20,20,21,55,56,110]# [-50,55]
##############################################################################
f_mn = fcts.fig1_multi(cutarn[0], cutarn[1], yla, nams, xla,12, hrzl)
sp = fcts.mkp(drc+'/testplts/'+Date,'/'+spec+'multi_norm'+cst+'_'+nm+'all')#curna)
f_mn.savefig(sp, dpi = 250)
#%% figure multi plot _ as measured
##############################################################################
ylab = '$dI/dV \ [G/G_{0}]$' # '$dG/dV \ [(G/G_{0})V^{-1}]$' # 'current [\u00b5A]' #'$G_{normalized}$' #  '$dI/dV \ [mG/G_{0}]$'
##############################################################################
f_mn = fcts.fig1_multi(cutar[0], cutar[1], ylab, legs=nams, \
                      xlab='Bias [mV]',bf=12)
sp = fcts.mkp(drc + '/testplts/'+Date,'/'+spec+'as_measured'+ cst+'_'+ curna) 
f_mn.savefig(sp, dpi = 250)
#%% stacked plots

print('plot spectra - stacked plots')

#%% figure stack normalized curves whole specs
##############################################################################
yla = 'G [a.u.]' # '$G_{normalized}$' #  '$dI/dV \ [mG/G_{0}]$'
##############################################################################
tit = ' '
fstn = fcts.fig_stack(cutarn,yla,nams,tit,xla,ydist = 0.7, hgt=8,bf =13)
sp = fcts.mkp(directory + '/testplts/'+Date+'/', '/stack_norm_' + curna) 
fstn.savefig(sp, dpi = 250)
#%% figure stack as measured
##############################################################################
ylab = '$dI/dV \ [G/G_{0}]$'
##############################################################################
fst = fcts.fig_stack(cutar,ylab,nams,tit,xla,ydist = 0.1,hgt = 6, bf=15)
sp = fcts.mkp(directory + '/testplts/'+Date+'/', '/stack_as_meas_' + curna) 
#fst.savefig(sp, dpi = 250)

#%%============================================================================


print('2d - histogramm and averaged spectra')




#%% 2d histogramms and averaging
##############################################################################                 
Tit = ['empty']
bins = 100
logy = 'no'
cmap = 'jet'
smo = 'middle' # smoothing grade 'rough', middle,'fine' or 'no'--> disable smoothing

############################################################################## 
#%% 2D hist Normalized
fhist = fcts.fig_2dHist(cutarn, Tit, [xla, yla], bins, cmap,logy)
sp = fcts.mkp(directory + '/testplts/'+Date+'/', '/2Dhist_norm' + curna+spec) 
fhist.savefig(sp, dpi = 250)
#%% 2D hist as measured
fhist = fcts.fig_2dHist(cutar, Tit, [xla, ylab], bins, cmap,logy)
sp = fcts.mkp(directory + '/testplts/'+Date+'/', '/2D_Hist' + curna+spec) 
fhist.savefig(sp, dpi = 250)

#=============================================================================
#%% averaged plot  normalized
fave = fcts.fig_average(cutarn, bins, xla, yla, logy, smo)
sp = fcts.mkp(directory+'/testplts/'+Date+'/','/average_spec_norm'+curna+spec) 
fave.savefig(sp, dpi = 250)
#%% averaged plot as measured
fave = fcts.fig_average(cutar, bins, xla, ylab, logy, smo)
sp = fcts.mkp(directory + '/testplts/'+Date+'/', '/average_spec' + curna+spec) 
fave.savefig(sp, dpi = 250)

#%% 



print('scatter plots')


#=============================================================================

# check if parameter is in al dicts (G_values for exapmple). returns names if parameter is not all dcts.
#-------------- >   parachecker(fids[0], 'G_values') <----------------



#%% set the coloring of the plots
##############################################################################
style = 'jet'
#############################################################################   
G_fac = fcts.propflt(fids, 'Gfac')


#%% scatter mplots 
"""
The names of the scatterplots are composed like this:
    
    date/savestring/xlongname _vs_ ylongname style .png
Date is automically generated
savestring is given with the variable 'svst'
The longnames are automatically chosen. In the 'multiple datapoints per sweep'
section if the variable 'ads' is a list of strings (not None), it takes that 
this list for the Longnames.

"""

#%%
##############################################################################
par= fcts.plot_arrays(['GI', 'ndr_r'], fids, nad)
save = False# True# 
svst = 'NDR' # savestring to further specify (if needed)
##############################################################################  
fx = fcts.fig_ps(par, svst, drc, Date, sav=save,ano=nams, fit='n',dpi=300)



#%% show a single property

xxx = fcts.propflt(nad, '  ')




#%% multi scatter plot 
##############################################################################
vya = fcts.plot_arrays(['swnr','GIr','swnr', 'GIp', 'swnr', 'GIn'], fids, nad)
paths = (directory, 'all')
xna, yna = 'Integrals', 'swnr'
##############################################################################

fmu = fcts.figmulti_sc(vya,si=40,tis=[],xln=xna,yln=yna,da=Date,pa=paths)


#%% Plot extrema Data


print('peak positions - extrema')   



#=============== Single  datapoint per curve plotting ========================
#%% Filter out 1 Datapoint per curve in certain Bias region ==================
##############################################################################
fra = (55,125) # filter range for peak filtering
pexy = [fcts.propflt(fids,'Gpe_x'), fcts.propflt(fids, 'Gpe_y')] # peak x,y
pra = 'Peak_R2_OC6_'+str(fra[0])+'-'+str(fra[1]) #peakrange
##############################################################################
didvfilt  = fcts.filt_from_peaks(pexy, fra, G_fac) 
dinam = np.array(nams)[didvfilt[3]] # names
dinam = dinam.tolist()
print('considered curves with one datapoint in range -fra-: \n\t', dinam)
#% load dictionary list again with the updated numbers (dinam):
fds = fcts.fbysn(fids, dinam) # filtered dictionaries (2nd time...)
nadp = fcts.new_dict2(fds, pexy, fra, G_fac, dinam) # new attribute dictionary peaks      
#%% peak  plots
##############################################################################                 
Gpexy = fcts.plot_arrays(['Mposn', 'Gpyr'], nadp, nadp) 
save =  False # True #
##############################################################################  
fx = fcts.fig_ps(Gpexy, pra, drc, Date, sav=save, ano=dinam, fit='n', dpi=300)


#=============================================================================
#%% print max min (just if that needs to be checked or known)
what = 'Gpxr'
print(min(fcts.propflt(nadp,what)), max(fcts.propflt(nadp,what)))

#%%




#========================= multiple datapoints per sweep =======================
#%% Plot scatter of multiple Datapoints per curve
##############################################################################
x_data = ['dGpa_x', 'dGda_x']
y_data = ['dGpa_y', 'dGda_y']
zdata = 'G0' # None             -> None if no bar is wanted
ads = None # adstring --> ['xlongname', 'ylongname']
##############################################################################
xyz_d = fcts.plot_conc(x_data, y_data, zdata, fids, ads)
ff = fcts.fig_ps(xyz_d, pra, drc, Date, sav=save, ano=dinam, fit='n', dpi=300)





#%%
############################################################################
bins = 50
save = False
############################################################################
xy_hist = fcts.plot_conc(x_data, y_data, None, fids, None)

(xy_hist, pra, drc, Date, sav=save, ano=dinam, fit='n', dpi=300)
fhi = fct.fig_ps(xy_hist, pra, drc, Date, sav=save, dpi=200, hist='yes',\
                 logy = 'no', bins=bins)




#%%


    




print('polycollection 3d plots')





#%% prepare json loaded data:
# prepare curves normalized and not
cut = [-80,80]

verts, vertst, hiG, loG, Vlim, cutar, ver = prepare_curves(evar, cut, None)
arrn = prepare_curves(evar, cut, 'ma')
vertsn, vertstn, hiGn, loGn, Vlim, cutarn, vern = arrn

#%% figure polycollection as measured
#curs = fnums
yla = '$G/G_{0}$ '
f_p = fcts.polyplot(verts, loG, hiG, Vlim, yla, fnums, 0, evar[2], azm = -60)
#f_p.savefig(directory + '/testplts/sweeps/3d_' + curna, dpi = 300)

#%% figure poly collection normalize
yla = '$G_{normalized}$'
fpn = fcts.polyplot(vertsn,loGn,hiGn,Vlim, yla, curs ,0, ficu[2], azm = -60)
#fpn.savefig(directory + '/testplts/sweeps/3d_norm_' + curna, dpi = 300)
#%% figure linecollection normalized
yla = '$G_{normalized}$'
minp = min(loGn)
flin = fcts.lineplot(vern, loGn, hiGn, Vlim, yla, curs, minp, evar[2], azm = -70)
#flin.savefig(directory + '/testplts/sweeps/3d_line_norm_' + curna, dpi = 300)

#%% figure linecollection as measured
yla = '$G/G_{0}$ '
fli = fcts.lineplot(ver, loG, hiG, Vlim, yla, curs, 0, evar[2], azm = -70)
#fli.savefig(directory + '/testplts/sweeps/3d_line_' + curna, dpi = 300)









