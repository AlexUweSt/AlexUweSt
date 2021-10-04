# -*- coding: utf-8 -*-
"""
Created on Wed Dec  9 15:52:12 2020

@author: strobe46
"""

import os
import fcts_OC_conc_ as fct

#%% init
# =============================================================================
# directory = 'C:/Users/Strobel/Nextcloud/Evaluation_data/' + \
#     'IETS_measurements_revised/2020_measurements/S3_2020/C60/'+ \
#         'E11/06_11/54/'
# =============================================================================
directory = 'C:/Users/Strobe46/ownCloud/Evaluation_data/' + \
    'IETS_measurements_revised/2020_measurements/S3_2020/C60/'+ \
        'E11/10_11/58/'
spath = directory + '/OC_conc_png'
# =============================================================================
# directory = 'C:/Users/strobe46/ownCloud/Evaluation_data/' + \
#     'IETS_measurements_revised/2020_measurements/S3_2020/C60/'+ \
#         'E11/11_11/59/'
# =============================================================================


#%% concatenate all curves, extract list
ocad, oca = fct.conc_OC(j,jj) # all curves [dict, list]
print(oca[4])
#%% figure oc all (no sweeps)
f_oc_a = fct.f_conc_oc(oca[0], oca[1], oca[2], psw = 'no', \
                       loca ='lower center', swpn = oca[3], cun = oca[5])
#%% figure all + sweeps (als) --> fig_conc_oc(concarr,tGswe,tbre,Gbre,psw,loca,swpn,cun, size=32):
f_oc_als = fct.f_conc_oc(oca[0], oca[1], oca[2], psw ='yes', \
     loca = 'lower center', swpn = oca[3], cun = oca[5], si=30, br = 'axv')
#%% save plot all with sweeps
f_oc_a.savefig(spath + '/__all_oc_cycles_nosweeps.png', dpi=450)
f_oc_als.savefig(spath + '/_all_oc_cycles_sweeps_axv.png', dpi=450)
#%%figure all temperature (special)
fat = fct.G_temperature(oca[0], oca[2][0], tlim = [4,37], lab = 'Gs')
#%%figure all temperature (static)
fat = fct.G_temperature(oca[0], oca[2][0], tlim = [4,37], lab = 'Gd')
#%%figure all temperature (lia --> differential)
fat = fct.G_temperature(oca[0], oca[2][0], tlim = [4,37], lab = '$G_{static}-G_{diff}$')
#%% save all temperature
fat.savefig(spath + '/__all_oc_cycles_Temp.png', dpi=310)
fat.savefig(spath + '/__all_oc_cycles_Temp_lia.png', dpi=310)
#%% decide what OCs to filter out
OCf = [2,2]
ocnf = str(min(OCf)) + '-' + str(max(OCf))
occd, occ = fct.choose_OCs(ocad, OCf)

#%%  figure choosen + sweeps (cs)
f_oc_cs = fct.f_conc_oc(occ[0], occ[1], occ[2],psw ='yes',loca ='lower left',\
                      swpn = occ[3], cun = occ[4], si=66)
#%% save chosen sweeps plot
f_oc_cs.savefig(spath + '/__oc_cycles_'+ ocnf +'.png', dpi=400)
#%% figure choosen no sweeps (c)
f_oc_c = fct.fig_conc_oc(occ[0], occ[1], occ[2], psw = 'no', loca = 'lower left',\
                     swpn = occ[3], cun = None)
#%% save fig chosen no sweeps
f_oc_c.savefig(spath + '/__oc_cycles_'+ ocnf +'_no_sweeps.png', dpi=450)

#%%figure chosen temperature (lia --> differential)
fat = fct.G_temperature(occ[0], occ[2][0], lab = 'Gd', tlim = [5,12])
#%% save all temperature
fat.savefig(spath + '/__chosen_oc_cycles_Temp.png', dpi=310)

#%%


#%% choose filter by time (smaller Window)
tlim = [160.04,180]#[97.9, 158]
tnf = str(tlim[0]) + '-' + str(tlim[1])
octd, occt = fct.choose_time(ocad, tlim)

#  figure choosen + sweeps + sweep numbers (csn)
f_oc_csn = fct.f_conc_oc(occt[0], occt[1], occt[2],psw = 'yes', \
         loca ='lower left',swpn = occt[3], cun=occt[4], si=105, tp = (0,-5))
#%%
f_oc_csn.savefig(spath + '/__oc_cycles_zoom_'+tnf+'.png', dpi=400)
#f_oc_csn.savefig(spath + '/__oc_cycles_zoom_'+tnf+'.svg')

#%% plot same as Temperature Plot
fat = fct.G_temperature(occt[0], occt[2][0], tlim = [4,12], lab = 'Gd')
fat.savefig(spath + '/__oc__zoom_Temperature_'+tnf+'.png', dpi=310)







