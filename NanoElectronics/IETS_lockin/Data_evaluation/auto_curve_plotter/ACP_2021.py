# -*- coding: utf-8 -*-
"""
Created on Thu Mar 18 15:12:51 2021

@author: strobe46

this programm is build upon the 'IETS auto curve plotter' programms
"""
# -*- coding: utf-8 -*-

import json
#from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
#from numpy import diff
import time
import os
import fcts_ACP_ as fct
import copy



# look for files

# %% ========================================================================
 ################# define variables and arrrays ##########################
root = '/home/filipk/Documents/Alex/Evaluation_data/IETS_measure/'
adpath = 'Ge_nanowire_IETS/04_03_2021/50nm/Wire_11/' #'2021/S1_2021/AuAu/Rp_200k/18_03/cold_sweeps_phase_new/'
##########################################################################
directory = os.path.join(root, adpath)
# filepaths: fullsweep, upsweep, downsweep, OC_file
fps_ = fct.filefinder(directory)
fps = fct.filesorter(fps_[1])
used_OC = []
dgdvlab = '$dG/dV \ [(G/G_{0})/V]$'
# %%============================================================================
# LockAmpl need for state 1 and 2 C
########### find files ############################
adic = {'LockAmpl': 0.004}  # Lock in AC signal amplitude
# adic['CAcorr'] = 1 # Current amplifier corection
# +- Factors to multiply IV GV and dGV curves (if amplifier is inverting or phase is shiftet -180 etc.)
adic['facs'] = (-1, -1, 1)
adic['Rpre'] = 2  # pre resistor for almost all measurements 2kOhm
adic['state'] = 5
adic['skiprows'] = 6
indx = 0
#################################################################

meas, fi_sw, sec_sw, header = fct.load_mainsweep(fps[0], adic, indx)

headers = fct.Header(header, adic)

full_arr = fct.calc_arrays(meas[0], headers[0], adi=adic)
fis_arr = fct.calc_arrays(fi_sw[0], headers[0], adi=adic)
sec_arr = fct.calc_arrays(sec_sw[0], headers[0], adi=adic)
dic_lst = [full_arr, fis_arr, sec_arr]
dlst_lst = fct.tolist(dic_lst)

hedict = fct.G_zero(full_arr, headers[0])
dilst = copy.deepcopy(dic_lst)
diclst, gylab, gfac = fct.G_scale(dilst, hedict)

fch = fct.figcheck(diclst[0], ['Bias_DMM', 'LIA1corr'], -100,
                       hedict['conductance']*0.84, gylab, 'amplitude ' + \
                           str(hedict['LockAmpl']))

dIkeys = fct.keymake(diclst[0], corr='yes')
num = fct.numeric(diclst[0], dIkeys, gfac,'middle')
Gzero = hedict['G_zero']*gfac
#%% numerics
fnum = fct.fig_deriv(diclst[0], dIkeys, num, gylab, Gzero,lim=15)

# %% loop 
allsweepdat =[]
for ind, path in enumerate(fps[0][0]):
    if os.path.getsize(path) > 41000:
        
        meas, fi_sw, sec_sw, header = fct.load_mainsweep(fps[0], adic, ind)
        headers = fct.Header(header, adic)  
        full_arr = fct.calc_arrays(meas[0], headers[0], adi=adic)
        fis_arr = fct.calc_arrays(fi_sw[0], headers[0], adi=adic)
        sec_arr = fct.calc_arrays(sec_sw[0], headers[0], adi=adic)
        dic_lst = [full_arr, fis_arr, sec_arr]
        dlst_lst = fct.tolist(dic_lst)
        hedict = fct.G_zero(full_arr, headers[0])
        dilst = copy.deepcopy(dic_lst)
        diclst, fct.gylab, gfac = fct.G_scale(dilst, hedict)
        
        allsweepdat.append([fct.sweepdata(ind, headers, fps[1][0][ind])])

        ############################# save data ##########################
        # define additional paths

        OCcyc = str(int(hedict['O/C_Cyles']))
        idx = str(ind)
        file = str(fps[1][0][ind]) + 'dia'
        OC_dir = directory + 'OC_' + OCcyc

        if os.path.exists(OC_dir) == False:
            os.mkdir(OC_dir)
        newdir = OC_dir + '/' + file
        if os.path.exists(newdir) == False:
            os.mkdir(newdir)
        mfp = newdir + '/Measurement_files'  # measurement file path
        if os.path.exists(mfp) == False:
            os.mkdir(mfp)
        paf = newdir + '/processed_arrays_files'
        if os.path.exists(paf) == False:
            os.mkdir(paf)
        plpa = newdir + '/additional_plots'  # plot path
        if os.path.exists(plpa) == False:
            os.mkdir(plpa)

        # save arrays:
        #np.savetxt(paf + '/V_IETS.txt' , V )
        np.savetxt(paf+'/fullsweep.txt',meas[1],header=','.join(meas[2]),\
                   fmt='%.6f')
        np.savetxt(paf + '/1st_sweep.txt', fi_sw[1],header =','.join(meas[2]),\
                   fmt='%.6f')
        np.savetxt(paf + '/2nd_sweep.txt', sec_sw[1],header=','.join(meas[2]),\
                   fmt='%.6f')
        if adic['state'] == 4:
            np.save(paf + '/time', time)
        np.savetxt(paf + '/Header.txt', headers[1], fmt='%s')

        nams = ['fullsweep', 'first_sweep', 'second_sweep']

        fct.save_dictionary([hedict], paf, ['header_dict'])
        fct.save_dictionary(dlst_lst, paf, nams)
        fct.copy_files([fps[0][0][ind], fps[0][1][ind], fps[0][2][ind]], mfp)

        # Single plots

        plt.close('all')

# ============================================================================================
        ###################### plots #################################################
        GOC = str(hedict['conductance'])
        if 'Temperature' in hedict.keys():

            Legend = GOC + ' G0 - '+str(hedict['Temperature']) + ' K'
        else:
            Legend = GOC + ' G0'

        fch.savefig(newdir + '/tstplt_dIdV_Goc')
        # plot stack:
        kys = fct.keymake(diclst)
        xlab = 'Bias [mV]'
        fst = fct.fig_stack(diclst[0], Legend, kys, xlab, gylab)  # fig stack
        fst.savefig(newdir + '/stack_black.png', dpi=300)

        # with first and second sweep

        fsud = fct.fig_up_down(diclst, gylab, kys)

        dIkeys = fct.keymake(diclst[0], corr='no')
        num = fct.numeric(diclst[0], dIkeys, gfac)
        Gzero = hedict['G_zero']*gfac
        fnum = fct.fig_deriv(diclst[0], dIkeys, num, gylab, Gzero, lim=20)
        fnum.savefig(plpa + '/numerical_comparison.png', dpi=300)


######################## fullsweep plots #####################################

        leg = str(hedict['Temperature']) + ' K'

        labels = ['Bias [mV]', 'Current [\u00b5A]']
        kes = ['Bias_DMM', 'I(DMM)']
        f1 = fct.fig_sing(diclst[0], kes, labels, leg)
        f1.savefig(newdir + '/IV_fullsweep.png', dpi=200)

        labels = ['Bias [mV]', gylab]
        kes = ['Bias_DMM', 'X_LIA1corr']
        f2 = fct.fig_sing(diclst[0], kes, labels, leg)

        if 'LIA1' or 'LIA_1' or 'LIA_1_DMM' in diclst[0].keys():
            labels = ['Bias [mV]', gylab]
            kes = ['Bias_DMM', 'LIA1corr']
            f3 = fct.fig_sing(diclst[0], kes, labels, leg)
            f3.savefig(newdir + '/LIA1_DMM_fullsweep.png', dpi=200)
            f2.savefig(plpa + '/LIA1_fullsweep.png', dpi=200)
        else:
            f2.savefig(newdir + '/LIA1_fullsweep.png', dpi=200)

        labels = ['Bias [mV]', dgdvlab]
        kes = ['Bias_DMM', 'LIA2_DMM']
        f4 = fct.fig_sing(diclst[0], kes, labels, leg)
        f4.savefig(newdir + '/LIA2_DMM_fullsweep.png', dpi=200)

        f5 = fct.fig_sing(diclst[0], ['Bias_DMM', 'X_LIA2'],
                      ['Bias [mV]', dgdvlab], leg)
        f5.savefig(plpa + '/LIA2_fullsweep.png', dpi=200)

        labels = ['Bias [mV]', 'IETS [1/V]']
        kes = ['Bias_DMM', 'XY']
        f6 = fct.fig_sing(diclst[0], kes, labels, leg)
        f6.savefig(newdir + '/IETS_fullsweep.png', dpi=200)

######################## allllsweep plots #####################################

        labels = ['Bias [mV]', 'Current [\u00b5A]']
        kes = ['Bias_DMM', 'I(DMM)']
        f7 = fct.fig_sing_updown(diclst, kes, labels)
        f7.savefig(plpa + '/IV_updown.png', dpi=150)

        labels = ['Bias [mV]', gylab]
        if 'LIA1' in diclst[0].keys():
            kes = ['Bias_DMM', 'LIA1corr']
        else:
            kes = ['Bias_DMM', 'X_LIA1corr']
        f8 = fct.fig_sing_updown(diclst, kes, labels)
        f8.savefig(plpa + '/dIdV_updown.png', dpi=150)

        labels = ['Bias [mV]', dgdvlab]
        kes = ['Bias_DMM', 'LIA2_DMM']
        f9 = fct.fig_sing_updown(diclst, kes, labels)
        f9.savefig(plpa + '/dGdV_updown.png', dpi=150)

        labels = ['Bias [mV]', dgdvlab]
        kes = ['Bias_DMM', 'XY']
        f9_1 = fct.fig_sing_updown(diclst, kes, labels)
        f9_1.savefig(plpa + '/IETS_updown.png', dpi=150)

        labels = ['Bias [mV]', gylab]
        kes = ['Bias_DMM', 'R_LIA1']
        f10 = fct.fig_sing_updown(diclst, kes, labels)
        f10.savefig(plpa + '/R_LIA1_updown.png', dpi=150)

        labels = ['Bias [mV]', 'Theta [degrees]']
        kes = ['Bias_DMM', 'thet_LIA1']
        f11 = fct.fig_scat_updown(diclst, kes, labels)
        f11.savefig(plpa + '/Thet_LIA1_updown.png', dpi=150)
        if 'R_LIA2' in diclst[0].keys():
            labels = ['Bias [mV]', dgdvlab]
            kes = ['Bias_DMM', 'R_LIA2']
            f12 = fct.fig_sing_updown(diclst, kes, labels)
            f12.savefig(plpa + '/R_LIA2_updown.png', dpi=80)

        labels = ['Bias [mV]', 'Theta [degrees]']
        kes = ['Bias_DMM', 'Th_LIA_2']
        f13 = fct.fig_scat_updown(diclst, kes, labels)
        f13.savefig(plpa + '/Thet_LIA2_updown.png', dpi=150)


######################### fullsweep zoom pics #################################

        limi = [-200, 200]

        labels = ['Bias [mV]', gylab]
        kes = ['Bias_DMM', 'X_LIA1corr']
        f14 = fct.fig_sing_lim(diclst[0], kes, labels, limi, leg)

        if 'LIA1' in diclst[0].keys():
            labels = ['Bias [mV]', gylab]
            kes = ['Bias_DMM', 'LIA1corr']
            f15 = fct.fig_sing_lim(diclst[0], kes, labels, limi, leg)
            f15.savefig(newdir + '/LIA1_DMM_fullsweep_zoom.png', dpi=200)
            f14.savefig(plpa + '/LIA1_fullsweep_zoom.png', dpi=200)
        else:
            f14.savefig(newdir + '/LIA1_fullsweep_zoom.png', dpi=200)

        labels = ['Bias [mV]', dgdvlab]
        kes = ['Bias_DMM', 'LIA2_DMM']
        f16 = fct.fig_sing_lim(diclst[0], kes, labels, limi, leg)
        f16.savefig(newdir + '/LIA2_DMM_zoom.png', dpi=200)


# %% OC curve filemaker

# G header: 'time', 'velocity', 'G', 'vdrop', 'Lia1', 'LIA2', 'G_LIA1', 'Sens', 'Gain', 'Acamp', 'temp'
# allsweepdat : allswehead = ['idx' , 'G0' , 'motpos' , 'data_points' , 'OC_cyc' , 'name']

Num = 0
for i in range(len(allsweepdat)):

    usd = allsweepdat[i][4]
    if usd not in used_OC:
        used_OC.extend([usd])

for count, dat in enumerate(fps[3]):
    if count in used_OC and os.path.getsize(dat) > 4500:

        #dct = load_oc(path)
        odc = fct.negvals(fct.load_oc(dat,adic))

        path_ = directory + 'OC_' + str(count) + '/'
        #np.savetxt(path_ + 'OC_file_' + str(count) + fps[1][3][dat], OC_file)

# ====  fig comparison
        fcompare = fct.fig_compare(odc, allsweepdat)
        fcompare.savefig(path_ + 'OC_file_' + str(count) + '-G_compared')

############################
        

        MposOC, Gval, LIAG = odc['MposOC'], odc['Gval'], odc['LIAG']
        swpos, tempe = ???? , odc['tempe']
        ti0 = 'Open close trace static conductance'
        fG0 = fct.G_colorjet_IV(MposOC, Gval, swpos, tempe,
                            [0], ti0, design='plasma')
        fG0.savefig(path_ + 'OC_file_' + str(count) + '-G_cmap')

        ti1 = 'Open close trace differential Conductance'
        fG1 = fct.G_colorjet_IV(MposOC, LIAG, swpos, tempe,
                            [0], ti1, design='plasma')
        fG1.savefig(path_ + 'OC_file_' + str(count) + '-LIA1_cmap')

    elif count not in used_OC:
        OCs_path = directory + 'OC_no_sweeps'
        if os.path.exists(OCs_path) == False:
            os.mkdir(OCs_path)
        OC_file = np.loadtxt(o[count], delimiter='\t', dtype=np.str)
        OC_file = np.char.replace(OC_file, ',', '.')
        OC_file = np.char.replace(OC_file, 'b', '').astype(np.float64)
        oc = np.transpose(OC_file)
        try:
            if len(oc[0]) > 40:
                timestamp = oc[0]
                MposOC = oc[1]/1E+7
                vel = oc[2]
                Gval_ = oc[3]
                Rval = oc[4]
                LIAG_ = oc[7]
                tempe = oc[11]

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

                plt.close('all')
                #path = directory + 'OC_' + str(count) + '/'
                np.savetxt(OCs_path + '/OC_file_' + dat, OC_file)

                ti2 = 'Open close trace static conductance'
                fG2 = G_colorjet(MposOC, Gval, tempe, ti2, design='plasma')
                fG2.savefig(OCs_path + '/' + dat + '-G_cmap.png')

                ti3 = 'Open close trace differential Conductance'
                fG3 = G_colorjet(MposOC, LIAG, tempe, ti3, design='plasma')
                fG3.savefig(OCs_path + '/' + dat + '-LIA_G_cmap.png')

                fig_compare(dct, allsweepdat)

                fig, ax = plt.subplots(figsize=(9.6, 4.8))
                Gp, = ax.scatter(MposOC, Gval, s=0.15, color='b', label='G')
                Gpd, = ax.scatter(MposOC, LIAG, s=0.15,
                                  color='m', label='diff G')
                if min(Gval) > min(LIAG):
                    plt.ylim(bottom=min(LIAG)-min(LIAG)/2,)
                else:
                    plt.ylim(bottom=min(Gval)-min(Gval)/2, )
                if max(Gval) > max(LIAG):
                    plt.ylim(top=max(Gval)+max(Gval)/2)
                else:
                    plt.ylim(top=max(LIAG)+max(LIAG)/2)
                ax.set_yscale('log')
                plt.grid(True, which='major', axis='y',
                         color='g', linewidth=0.2)
                ax.set_xlabel('Motorposition [a.u.]', fontsize=20)
                ax.set_ylabel('conductance $[G/G_0]$', fontsize=20)
                ax.tick_params(axis='both', labelsize=18.0)
                ax2 = ax.twinx()
                ax2.scatter(MposOC, vel, s=0.15, color='k', alpha=0.15)
                plt.legend(
                    [Gp, Gpd], ['static G', 'differential G'], loc='best')
                ax2.tick_params(axis='y', labelsize=16.0)
                ax.set_ylabel('Motorspeed', fontsize=16)
                if abs(MposOC[0]) > abs(MposOC[len(MposOC)-1]):
                    ax.invert_xaxis()
                plt.tight_layout()
                # plt.show()
                plt.savefig(OCs_path + '/' + dat + '-G_plot.png')

        except TypeError:
            pass

#np.savetxt(directory + '/Allsweepdata.txt', allsweepdat, fmt='%s')


