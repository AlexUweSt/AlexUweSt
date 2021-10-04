# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 09:45:46 2020

@author: strobe46
"""

import scipy.constants as spc
import matplotlib.pyplot as plt
import numpy as np

def fig1_single(T,E,unit,temp):
    """ signal 'iets' or 'dGdV', I or didv """
    if temp == 'K':
        xlabel = 'Temperature [K]'
    elif temp == 'mK':
        xlabel = 'Temperature [mK]'
    f1 = plt.figure(1,figsize=[6.4, 4.8])
    ax1 = f1.add_subplot(111)
    #ax4 = ax3.twinx()
    f1.tight_layout()
    if unit == 'eV':
        ax1.plot(T,E, c='K')
        ax1.set_ylabel('Energy [eV]', fontsize = 16)
    elif unit == 'Joule':
        ax1.set_ylabel('Energy [Joule]', fontsize = 16)
        ax1.plot(T,E, c='K')
    elif unit == 'meV':
        ax1.set_ylabel('Energy [meV]', fontsize = 16)
        ax1.plot(T,E, c='K')
    elif unit == 'mJoule':
        ax1.set_ylabel('Energy [mJole]', fontsize = 16)
        ax1.plot(T,E, c='K')
    ax1.set_xlabel(xlabel, fontsize = 16)
    ax1.tick_params(axis='both' , labelsize = 12.0)
    #ax1.legend(['T = ' + str(temp)+ ' K'])
    f1.tight_layout()
    return f1

def mev_cm_converse(wvnr):
    meV = wvnr*spc.c*100*spc.h/(spc.e/1E3)
    return meV

def cm_meV_converse(meV):
    wvnr = meV/(spc.c*100*(spc.h/spc.e))/1E3
    return wvnr

def meV_THZ_converse(meV):
    THz = meV*spc.e/(spc.h*1E15)
    return THz
#%% constants
kb = spc.Boltzmann
e = spc.eV
me = spc.electron_mass
h = spc.h
hb = spc.hbar
fermi = spc.fermi
EeV = spc.eV

#%% Temperature - energy conversation
# arrays

Thigh = np.arange(1, 300, 2)
Tlow = np.arange(1, 25, 5E-1)

TlowmK = Tlow*1000

#%% calc Energies

Ejh = kb*Thigh

Ejl = kb*Tlow

Emjh = kb*Thigh*1000

Emjl = kb*Tlow*1000

EeVh = kb * Thigh / e

EeVl = kb * Tlow / e

EmeVh = (kb * Thigh / e)*1000

EmeVl = (kb * Tlow / e) * 1000

#%% Plot

f1 = fig1_single(Thigh,Ejh,'Joule','K')

f2 = fig1_single(Thigh,EeVh,'eV','K')

f3 = fig1_single(Thigh,EmeVh,'meV','K')

f4  = fig1_single(TlowmK,EmeVl,'meV','mK')

f5 = fig1_single(Tlow,EmeVl,'meV','K')

#%% Energy unit conversation




#%% Energy wavenumber conversation
###############################################################################
wavenr = 300#1492#496
###############################################################################
meV = mev_cm_converse(wavenr)
print('\n'+str(meV)+ ' meV')
#%%  wavenumber Energy conversation
###############################################################################
energ = 37.2  # in meV
###############################################################################
shift = cm_meV_converse(energ)
print('\n'+str(shift)+ ' cm^-1')

#%% energy THz converse
############################################################################
energ = 37.2
##########################################################################
THz = meV_THZ_converse(energ)
print('\n'+str(THz)+ ' THz')
#%% Energy





#%% siemens G0 conversation
###############################################################################
siemens = 1.5E-5
###############################################################################
G0 = spc.physical_constants['conductance quantum'][0]
gg0 = siemens/G0
print('\n' + str(round(gg0,5)) + ' G/G0')


#%% Fermi conversations
##############################################################################
efm = 0.055 # effective mass factor
Energ = 0.3 # in eV
fac = 1E6
###############################################################################
#Wave Vector Fermi
meff = efm*me
EJ = Energ*EeV
kF = ((2*meff*EJ)/hb**2)**0.5 # Fermi vector

# Fermi velocity

vF = (hb*kF)/meff
vF_ = ((2*EJ)/meff)**0.5

print(vF/fac, vF_/fac)

#%% Fermi Energy from velocity
###############################################################################
vFgiven = 1.4E6
efm = 0.055 # effective mass factor
###############################################################################
EFcalc = ((vFgiven**2*efm*me)/2)/EeV

print(EFcalc)
#%% lenght Energy Conversation
###############################################################################
efm = 0.051 # effective mass factor
Energ = 0.3 # Fermi energy in eV
fac = 1E6 # randmom factor
Ec = 0.25 # thouless energy scale in ev
###############################################################################
#Wave Vector Fermi
meff = efm*me
EJ = Energ*EeV
kF = ((2*meff*EJ)/hb**2)**0.5 # Fermi vector
# thouless energy in Joule
EcJ = Ec*EeV
#Fermi velocity
vF = (hb*kF)/meff
vF_ = ((2*EJ)/meff)**0.5
#print(vF/fac, vF_/fac)

# l calculation

lth = (hb*((2*EJ)/meff)**0.5)/(EcJ)

lc = hb*vF/EcJ
print('\n')
print(lth*1E9,lc*1E9)
#print(lth,lc)



















