#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Library to capture rate calculations
author : arturonunez25@gmail.com
"""
import numpy as np
import os
from numpy import exp, sqrt
#import BORZfits as dm
import fdulib as vdf
import math
import get_nuclear_info
import ConfigParser
from wkbl.inputing import dm_input

######## getting input ########
D = dm_input()
D.parse_input(os.getcwd()+"/input.dm")

#######################################
################ SUN ##################

sun = get_nuclear_info.Sun

def say_Sun():
    return sun

#######################################

def get_var():
    """
    """
    v_Sun = float(raw_input("Sun velocity around the GC\nv_Sun = "))
    v_esc = float(raw_input("Escape Velocity of the WIMPs in the HAlO\nv_esc = "))
    v_0 = float(raw_input("Dispersion speed\nv_0 = "))
    rho = float(raw_input("Local DM density\nrho ="))
    os.system('cls' if os.name == 'nt' else 'clear')
    return v_Sun, v_esc, v_0, rho


def get_variables():
    """
    prints current variables in use
    """
    print "Solar velocity arround GC     v_Sun = ", D.v_Sun," km / s"
    print "Wimp density                    rho = ", D.rho," GeV / km^3 "
    print "Wimp-nucleon cross section    sigma = ", D.sigma_sd, "km^2 SD Wimp-proton    \nffrom Ellis,Olive and Savage (2008) "
    print " proton mass                    m_p = ", D.m_p, " GeV for SD is Hydorgen mass"
    print "Nuclear radius(depends on m_i)  r_i = ", D.r_h, " ^km for hydrogen"
    print "Escape velocity               u_esc = ", D.v_esc, "km / s"
    print "Dispertion velocity             v_0 = ", D.v_0, "km /s"
    return 0

def movingWA(lengh,x_array, array):
    """
    moving window average
    for two array x and y
    """
    x = [0 for i in range(lengh)]
    y = [0 for i in range(lengh)]
    x_final = []
    y_final = []

    for i in range(0,lengh):
            if math.isnan(array[i])== True:
                    continue
            x[i%lengh] = x_array[i]
            y[i%lengh] = array[i]
    for i in range(lengh+1 , len(array)):
            if math.isnan(array[i])== True:
                    continue
            x_final.append(sum(x)/lengh)
            y_final.append(sum(y)/lengh)
            x[i%lengh] = x_array[i]
            y[i%lengh] = array[i]
    return x_final, y_final

def normalize(array):
    """
    normalizes an array to a const value
    inside the array
    """
    const = array[1]
    nuarray = []
    for i in range(len(array)):
        nuarray.append(array[i] / const)
    return nuarray

def normalize_max(array):
    """ 
    normalizes an array to a const value
    inside the array
    """
    const = np.amax(array)
    nuarray = []
    for i in range(len(array)):
        nuarray.append(array[i] / const)
    return nuarray

def exp_ax2(x , a):
    """
    gaussian if a is negative
    """
    return exp(a * (x**2))

def compare(f1,m,fit):
    """

    :param f1: simulation to be used
    :param m: mass of the wimp
    :param fit: fit fuction tu be used
    :return:
    """
    val = quad(product, 0., 1100., args=(dm.fdu['maxwellian'], capLib.caprate_GOU,
                                         m, f1['maxwellian'], 270., 1080., sun[0]))
    val1 = quad(product, 0., 1100., args=(dm.fdu[fit], capLib.caprate_GOU,
                                          m, f1[fit], 270., 1080., sun[0]))
    if val[0]==0.:
        return float('nan')
    else:
        return val1[0] / val[0]


def get_arg(num):
    if num == 1:
        fun = dm.e1
    elif num == 2:
        fun = dm.e2
    elif num == 3:
        fun = dm.e3
    elif num == 4:
        fun = dm.e4
    elif num == 5:
        fun = dm.e5
    elif num == 6:
        fun = dm.e6
    elif num == 7:
        fun = dm.e7
    elif num == 8:
        fun = dm.e8
    elif num == 9:
        fun = dm.e9
    elif num == 10:
        fun = dm.e10
    elif num == 11:
        fun = dm.e11
    elif num == 12:
        fun = dm.e12
    else:
        print "ERROR wrong sim number"
        sys.exit()

    return fun




def caprate_GOU(u, m , rms='s'):
    """
    this is the gould expression
    input
    u = velocity
    m = WIMP mass
    u_e = escape velocity from the sun
    v_r = dispertion velocity
    """
    if rms == 's':
            v_r = D.v_shm
    elif rms == 'l':
            v_r = D.v_lin
    elif rms == 'm':
            v_r = D.v_mao
    else:
        try:
            v_r = np.float(rms)
        except:
            print "ERROR: rms velocity not properly defined"
            print "           rms = ", rms, "\n\n"
    u_e = D.v_esc #+ v_Sun### check
    nuc = sun[0]
    # information from the nucleus
    g =float(nuc[3])
    m_i = float(nuc[2])
    ru =(m_i * m)/(m_i + m)  # WIMP nucleous reduced mass
    ru_p = (D.m_p * m)/(D.m_p + m)  # WIMP proton reduced mass
    # SI 0-momentum transfer cross section (CHOI)
    sigma_i = (D.sigma_sd * (float(nuc[0])**2) * ru) / ru_p
    # mu factor eq 2.10 GOULD 87
    mu = m_i / m
    mu_p = (mu + 1.) / 2.
    mu_m = (mu - 1.) / 2.
    x = u / v_r
    # factor A^2 eq 2.22 GOULD 87
    # 21 is a empiric number might be wrong
    k2 =(21**2) * g * (mu / (mu_m**2))
    # original form is
    # k2 =(3./2.) * ((u_e**2) / (v_r**2)) * (mu / (mu_m**2))
    # small a factor and small b factor eq A7 GOULD 87

    a = 1e13 * (2./9.) * m * m_i * ((v_r * D.r_h)**2)
    b = (mu * a) / (mu_p**2)
    # final expression eq A6 GOULD 87
    # constant term
    const = 1e36 * (D.rho * sigma_i) / (m * b)
    # u dependent terms
    c_1 = exp(-a * m * (x**2))
    c_2_1 = exp(-k2 * (a-b)) # aA2 depends on u
    c_2_2 = exp(-b * m * (x**2))
    # --flag--
    # print "a =", a, " k =",k2, "m =", m, " m_p =", m_p," const", const
    return const * (c_1 - (c_2_1 * c_2_2))

def caprate_GOUSI(u, m, index):
    """
    this is the gould expression
    input
    u = velocity
    m = WIMP mass
    u_e = escape velocity from the sun
    v_r = dispertion velocity
    """
    nuc = sun[index]
    # information from the nucleus
    u_e = D.v_esc
    v_r = D.v_shm
    g = float(nuc[3])
    m_i = float(nuc[2])
    ru =  (m_i * m)/(m_i + m) # WIMP nucleous reduced mass
    ru_p = (D.m_p * m)/(D.m_p + m) # WIMP proton reduced mass
    # SI 0-momentum transfer cross section (CHOI)
    sigma_i =  (D.sigma_si * (float(nuc[0])**2) * ru) / (ru_p)
    # mu factor eq 2.10 GOULD 87
    mu = m_i / m
    mu_p = (mu + 1.) / 2.
    mu_m = (mu - 1.) / 2.
    x = u / v_r
    # factor A^2 eq 2.22 GOULD 87
    k2 =((21)**2) * (g) * (mu / (mu_m**2))
    #k2 =(3./2.) * ((u_e**2) / (v_r**2)) * (mu / (mu_m**2))
    # small a factor and small b factor eq A7 GOULD 87
    a = (1e13) * (2./9.) * m * m_i * ((v_r * D.r_h)**2)
    b = (mu * a) / (mu_p**2)
    # final expression eq A6 GOULD 87
    # constant term
    const = 1e36 * float(nuc[1]) * (D.rho * sigma_i) /(m * b)
    # u dependent terms
    c_1 = exp(-a * m * (x**2))
    c_2_1 = exp(-k2 * (a-b)) # aA2 depends on u
    c_2_2 = exp(-b * m * (x**2))
    # --flag--
    # print "a =", a, " k =",k2, "m =", m, " m_p =", m_p," const", const
    result = const * (c_1 - (c_2_1 * c_2_2))
    if result > 0:
        return result
    else:
        return 0






