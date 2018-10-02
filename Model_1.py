#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 20:11:47 2018

@author: shilpakancharla
"""

import math
import matplotlib.pyplot as plt
from scipy.special import gamma

#Functions for converting units - convert to cgs units

# -*- coding: utf-8 -*-
"""
For information and usage see README, or http://pypi.python.org/pypi/numericalunits
"""
#Copyright (C) 2012-2017 Steven Byrnes
#
#Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
#
#The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
#
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#From the numericalunits package: https://pypi.org/project/numericalunits/

#Set all variables to help introspection libraries
m = kg = s = C = K = 0.

cm = mm = um = nm = pm = fm = km = angstrom = lightyear = \
    astro_unit = pc = kpc = Mpc = Gpc = inch = foot = mile = thou = 0.
    
ms = us = ns = ps = fs = minute = hour = day = week = year = 0.
    
J = mJ = uJ = nJ = pJ = fJ = kJ = MJ = GJ = erg = eV = meV = keV = MeV = GeV = \
    TeV = btu = smallcal = kcal = Wh = kWh = 0.
    
def reset_units(seed=None):
    
    #Set all units to new self-consistent, floating-point values. 
    
    """
    reset_units() --> units are randomized. This is the suggested use. Run this 
    before your calculation, display the final answer, then re-run this, then re-display
    the final answer. If you get the same answers both times, then your calculations
    are almost guaranteed to be free of dimensional-analysis-violating errors.
    This method is run automatically the first time this module is imported.
    
    reset_unit('SI') --> Set units so that all values are given in standard SI
    units (meters-kilograms-seconds) by default. In this mode, there is no way to test
    for dimensional-analysis-violating errors.
    
    reset_units(x) --> If you pass any other argument x, it's used as the seed
    for the random number generator.
    
    """
    
    import random
    
    global m, kg, s, C, K
    
    if seed == 'SI':
        m = 1.
        kg = 1.
        s = 1.
        C = 1.
        K = 1.
    else: 
        prior_random_state = random.getstate()
        
        if seed is Non: 
            random.seed()
        else:
            random.seed(seed)
        
        m = 10. ** random.uniform(-1, 1) #meter
        kg = 10. ** random.uniform(-1, 1) #kilogram
        s = 10. ** random.uniform(-1, 1) #second
        C = 10. ** random.uniform(-1, 1) #coulombs
        K = 10. ** random.uniform(-1, 1) #kelvins
        
        #Leave the random generator like found, in case someone else is using it.
        random.setstate(prior_random_state)
    
    set_derived_units_and_constants()
    return

def set_derived_units_and_constants():
    
    """
    Assuming that the base units (m, kg, s, C, K) have already been set as 
    floating point values, this function sets all other units and constants to the
    appropriate, self-consistent values.
    
    """
    
    #Length
    global cm, mm, km, lightyear, astro_unit, pc, kpc, Mpc, Gpc, inch, foot, mile
        
    cm = 1e-2 * m
    mm = 1e-3 * m
    km = 1e3 * m
    lightyear = 9460730472580800. * m
    astro_unit = 149597870700. * m #astronomical unit
    pc = (648000./math.pi) * astro_unit #parsec
    kpc = 1e3 * pc
    Mpc = 1e6 * pc
    Gpc = 1e9 * pc
    inch = 2.54 * cm
    foot = 12. * inch
    mile = 5280. * foot
    
    #Time
    global minute, hour, day, week, year
    
    minute = 60. * s
    hour = 60. * minute
    day = 24. * hour #solar day
    week = 7. * day
    year = 365.256363004 * day #sidereal year
    
    #Energy
    global J, erg, eV, keV, MeV, GeV, Tev
    
    J = (kg * m**2) / s**2
    erg = 1e-7 * J
    eV = 1.6021766208e-19 * J
    keV = 1e3 * eV
    Mev = 1e6 * eV
    GeV = 1e9 * eV
    TeV = 1e12 * eV

reset_units('SI')
set_derived_units_and_constants()
    
#Distance between supernova and detector on Earth
#Using 51.4 kpc for SN1987A
dist = 10 * kpc

#Enumeration of matter type and flavor of neutrino

matterType = ['nu', 'nubar']
flavor = ['e', 'mu', 'tau']

#Luminosity
L_array = [[1.6e52 * (erg / s), 1.6e52 * (erg / s)],
            [1.6e52 * (erg / s), 1.6e52 * (erg / s)],
            [1.6e52 * (erg / s), 1.6e52 * (erg / s)]]
        
#Mean energy, 
meanE_array = [[0.015 * GeV, 0.015 * GeV],
                [0.015 * GeV, 0.015 * GeV],
                [0.015 * GeV, 0.015 * GeV]]
        
#Alpha (pinch parameters)

alpha_array = [[3, 4.31],
                [6, 6],
                [6, 6]]
                
#Creating a dictionary that is enumerated

alpha = {}
L = {}
meanE = {}
for x, pair_item in enumerate(flavor):
    for y, item in enumerate(matterType):
        alpha[(item, pair_item)] = alpha_array[x][y]
        L[(item, pair_item)] = L_array[x][y]
        meanE[(item, pair_item)] = meanE_array[x][y]

#Float range

def frange(start, stop, step = 0.0002):
    while start < stop:
        yield start
        start += step

def nuE(E):
    
    #Ratio of mean energy and luminosity
    ratio_e = (L[('nu', 'e')] / meanE[('nu', 'e')])
            
    #Numerator
    num_e = (alpha[('nu', 'e')] + 1) ** (alpha[('nu', 'e')] + 1)
            
    #Denominator
    den_e = (meanE[('nu', 'e')]) * gamma(alpha[('nu', 'e')] + 1)
            
    #Fraction
    frac_e = (num_e / den_e)
      
    flux_e = ((ratio_e) * (frac_e) * ((E / meanE[('nu', 'e')]) ** (alpha[('nu', 'e')])) * math.exp(-((alpha[('nu', 'e')] + 1) * E) / meanE[('nu', 'e')])) / (4*math.pi*dist**2)
    return flux_e
    
def nubarE(E):
    
    #Ratio of mean energy and luminosity
    ratio_anti_e = (L[('nubar', 'e')] / meanE[('nubar', 'e')])
            
    #Numerator
    num_anti_e = (alpha[('nubar', 'e')] + 1) ** (alpha[('nubar', 'e')] + 1)  
            
    #Denominator
    den_anti_e = (meanE[('nubar', 'e')]) * gamma(alpha[('nubar', 'e')] + 1)
            
    #Fraction
    frac_anti_e = (num_anti_e / den_anti_e)

    flux_anti_e = ((ratio_anti_e) * (frac_anti_e) * ((E / meanE[('nubar', 'e')]) ** (alpha[('nubar', 'e')])) * math.exp(-((alpha[('nubar', 'e')] + 1) * E) / meanE[('nubar', 'e')])) / (4*math.pi*dist**2)
    return flux_anti_e
    
def nuMu(E):
    
    #Ratio of mean energy and luminosity
    ratio_mu = (L[('nu', 'mu')] / meanE[('nu', 'mu')])
            
    #Numerator
    num_mu = (alpha[('nu', 'mu')] + 1) ** (alpha[('nu', 'mu')] + 1)
            
    #Denominator
    den_mu = (meanE[('nu', 'mu')]) * gamma(alpha[('nu', 'mu')] + 1)
            
    #Fraction
    frac_mu = (num_mu / den_mu)

    flux_mu = ((ratio_mu) * (frac_mu) * ((E / meanE[('nu', 'mu')]) ** (alpha[('nu', 'mu')])) * math.exp(-((alpha[('nu', 'mu')] + 1) * E) / meanE[('nu', 'mu')])) / (4*math.pi*dist**2)
    return flux_mu
    
def nubarMu(E):
    
    #Ratio of mean energy and luminosity
    ratio_anti_mu = (L[('nubar', 'mu')] / meanE[('nubar', 'mu')])
            
    #Numerator
    num_anti_mu = (alpha[('nubar', 'mu')] + 1) ** (alpha[('nubar', 'mu')] + 1)
            
    #Denominator
    den_anti_mu = (meanE[('nubar', 'mu')]) * gamma(alpha[('nubar', 'mu')] + 1)
            
    #Fraction
    frac_anti_mu = (num_anti_mu / den_anti_mu)
      
    flux_anti_mu = ((ratio_anti_mu) * (frac_anti_mu) * ((E / meanE[('nubar', 'mu')]) ** (alpha[('nubar', 'mu')])) * math.exp(-((alpha[('nubar', 'mu')] + 1) * E) / meanE[('nubar', 'mu')])) / (4*math.pi*dist**2)
    return flux_anti_mu
    
def nuTau(E):
    
    #Ratio of mean energy and luminosity
    ratio_tau = (L[('nu', 'tau')] / meanE[('nu', 'tau')])
            
    #Numerator
    num_tau = (alpha[('nu', 'tau')] + 1) ** (alpha[('nu', 'tau')] + 1)  
            
    #Denominator
    den_tau = (meanE[('nu', 'tau')]) * gamma(alpha[('nu', 'tau')] + 1)
            
    #Fraction
    frac_tau = (num_tau / den_tau)
    
    flux_tau = ((ratio_tau) * (frac_tau) * ((E / meanE[('nu', 'tau')]) ** (alpha[('nu', 'tau')])) * math.exp(-((alpha[('nu', 'tau')] + 1) * E) / meanE[('nu', 'tau')])) / (4*math.pi*dist**2)
    return flux_tau
    
def nubarTau(E):
    
    #Ratio of mean energy and luminosity
    ratio_anti_tau = (L[('nubar', 'tau')] / meanE[('nubar', 'tau')])
    
    #Numerator
    num_anti_tau = (alpha[('nubar', 'tau')] + 1) ** (alpha[('nubar', 'tau')] + 1)
            
    #Denominator
    den_anti_tau = (meanE[('nubar', 'tau')]) * gamma(alpha[('nubar', 'tau')] + 1)
            
    #Fraction
    frac_anti_tau = (num_anti_tau / den_anti_tau)
    
    flux_anti_tau = ((ratio_anti_tau) * (frac_anti_tau) * ((E / meanE[('nubar', 'tau')]) ** (alpha[('nubar', 'tau')])) * math.exp(-((alpha[('nubar', 'tau')] + 1) * E) / meanE[('nubar', 'tau')])) / (4*math.pi*dist**2)
    return flux_anti_tau
    
#File-output for energies of neutrinos (tab delimited) - will output from for loop
#Will also output graphs of neutrino spectra
#Can change start, stop, and step later - first testing to see if for loop works and writes a text file

for t in range(0, 20, 1):
    #Functions go here
    file = open("data/data_%d.txt" % t, "w+")
    for E in frange(0, 0.1002, 0.0002):
        flux_e = (nuE(E * GeV) * (s * cm**2) * 0.0002 * GeV) 
        flux_anti_e = (nubarE(E * GeV) * (s * cm**2) * 0.0002 * GeV)
        flux_mu = (nuMu(E * GeV) * (s * cm**2) * 0.0002 * GeV)
        flux_anti_mu = (nubarMu(E * GeV) * (s * cm**2) * 0.0002 * GeV)
        flux_tau = (nuTau(E * GeV) * (s * cm**2) * 0.0002 * GeV)
        flux_anti_tau = (nubarTau(E * GeV) * (s * cm**2) * 0.0002 * GeV)
        #print(str(E)+ '\t' + str(flux_e) + '\t' + str(flux_anti_e) + '\t' + str(flux_mu) + '\t' + str(flux_anti_mu) + '\t' + str(flux_tau) + '\t' + str(flux_anti_tau))
        stringData = str(E)+ '\t' + str(flux_e) + '\t' + str(flux_anti_e) + '\t' + str(flux_mu) + '\t' + str(flux_anti_mu) + '\t' + str(flux_tau) + '\t' + str(flux_anti_tau)
        file = open("data/data_%d.txt" % t, "a") 
        file.write(stringData + '\n')
        file.close()
    file = open("data/data_%d.txt" %t)
    lines = file.readlines()
    E = []
    flux_e = []
    flux_anti_e = []
    flux_mu = []
    flux_anti_mu = []
    flux_tau = []
    flux_anti_tau = []
    for line in lines:
        E.append(line.split()[0])
    for line in lines:
        flux_e.append(line.split()[1])
    for line in lines:
        flux_anti_e.append(line.split()[2])   
    for line in lines:
        flux_mu.append(line.split()[3])
    for line in lines:
        flux_anti_mu.append(line.split()[4])
    for line in lines:
        flux_tau.append(line.split()[5])
    for line in lines:
        flux_anti_tau.append(line.split()[6])
    file.close()
    plt.figure()
    plt.plot(E, flux_e, 'pink', lw = 2, label="Electron neutrino")
    plt.plot(E, flux_anti_e, 'mediumpurple', lw = 2, label="Electron antineutrino")
    plt.plot(E, flux_mu, 'lightskyblue', lw = 2, label="Mu neutrino")
    plt.plot(E, flux_anti_mu, 'silver', lw = 2, label="Mu antineutrino")
    plt.plot(E, flux_tau, 'gold', lw = 2, label="Tau neutrino")
    plt.plot(E, flux_anti_tau, 'red', lw = 2, label="Tau antineutrino")
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.title('Neutrino Spectra_%d' %t)
    plt.xlabel("Energy")
    plt.ylabel("Flux")
    plt.grid()
    plt.savefig('spectra/Neutrino Spectra_%d.png' %t, bbox_inches='tight')