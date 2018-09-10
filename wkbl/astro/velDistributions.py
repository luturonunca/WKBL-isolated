#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

Library of different velocity ditribution functions
f(u)/u
for dark halos
author : arturonunez25@gmail.com
"""
import sys
import math
import glob
import cmath
import emcee
import subprocess
import numpy as np
from unsio import *
import scipy.special as sp
from numpy import exp, sqrt
import scipy.integrate as integrate
from wkbl.inputing import dm_input
quad = integrate.quad

############################################################
################## getting input ###########################

D = dm_input()
D.parse_input(os.getcwd()+"/input.dm")

############################################################
############################################################


############################################################################
############################################################################
################# Velocity Distributio Functions ###########################
############################################################################
###############  in the galatic halo frame of reference  ###################
#normalization factors
# shm
f = lambda v: (v**2)*(1/(D.v_shm**3))*exp((-3.*(v**2))/(2. * D.v_shm**2 ))  
n_shm = 1./quad(f, 0., D.v_esc)[0]
n_shm1 = 1./quad(f, 0., D.v_esc2)[0]
n_shm2 = 1./quad(f, 0., D.v_esc3)[0]

p= 0.78 # mao fitting parameter
# mao
def f(v,v_e):
        return (v**2) * exp(-v/D.v_mao) * (((v_e**2) - (v**2))**p)
n_mao = 1/quad(f, 0., D.v_esc, args=(D.v_esc))[0]
n_mao1 = 1/quad(f, 0., D.v_esc2, args=(D.v_esc2))[0]
n_mao2 = 1/quad(f, 0., D.v_esc3, args=(D.v_esc3))[0]
# lin
q = 0.77
p =  q/(1-q)
def f (v):
    if v < 555:
            return (v**2) * ((1-((1-q) * (v**2 / D.v_lin**2)))**p)
    else:
            return 0
n_lin = 1/quad(f, 0., D.v_esc)[0]
n_lin1 = 1/quad(f, 0., D.v_esc1)[0]
n_lin2 = 1/quad(f, 0., D.v_esc2)[0]
###########################################################################


def maxw(u, v_S, v_r, u_esc):
    """
    Maxwellian VDF
    input: u variable velocity
           m WIMP mass
           v_S velocity of the Sun around the GC
           v_r dispertion velocity of the halo

    """
    u_esc = D.v_esc
    v_r = D.v_shm
    if u > u_esc:
            return 0
    else:
        a_f = (1./(v_S * v_r * sqrt(np.pi))) * exp(-(v_S**2) / (v_r**2)) # constant term
        b_f = u * exp(-(u**2) / (v_r**2))  #first term depending on u
        c_f = lambda u: exp(2. * u * v_S/(v_r**2)) # second (and third) term depending on u
        f = a_f * b_f * (c_f(u) - c_f(-u)) # maxwellian VDF in the Sun R. F
        return f



def shm_gal(v, transform=False, v_e=D.v_esc, n=n_shm):
        """
        maxwellian in the halo frame
        when esc diferent from zero 
        v_esc will be the provided
        value for esc
        """
        v_r = D.v_shm #dispersion speed
        if v > v_e:
                return 1e-30
        else:
                t = (1/(v_r**3)) * exp((-3. * (v**2)) / (2. * v_r**2 ))
                if transform == False:
                        return (n * (v**2) * t)
                if transform == True:
                        return  n * t

def mao_gal(v, transform=False, p=0.78, v_e=D.v_esc, n=n_mao):
        """
        Mao et all empiric expression
        in the halo frame of reference
        """
        v_r = D.v_mao #dispertion speed 
        if v > v_e:
                return 1e-30
        else:
                t1 = exp(-v/v_r)
                t2 = ((v_e**2) - (v**2))**p
                if transform == False:
                        return n * (v**2) * t1 * t2
                if transform == True:
                        return n * t1 * t2

def lin_gal(v, transform=False, v_e=v_esc, n=n_lin):
        """
        tsallis expression in the
        halo frame of reference
        from ling et al.
        """
        v_r = v_lin
        if v > v_e:
                return 1e-30
        else:
                t = (1-((1-q) * (v**2 / v_r**2)))
                if transform == False:
                        return n * (v**2) * t**p
                if transform == True:
                        return n * t**p

############################################################################
############ transforming it to the Sun frame of references ################

def transform(x, u, f, es, sun):
        """
        this is the galilean transformations that corresponds to the
        expression 2.2  in the CHOI paper; the part inside the integral
        """
        v_S = D.v_Sun
        if sun !=0:
                v_S = sun
        v = sqrt((u**2) + (v_S**2) + (2 * u * v_S * x)) 
        return f(v,transform = True, esc=es)

def gal_to_sun(u, f, esc_sun=[D.v_esc, D.v_Sun]): 
        """ 
        This is the transformation of a generic VDF in the galactic halo frame of 
        reference to the Sun f.o.r where
        v_S: speed of the Sun
        v_r: dispersion speed of the halo model
        v_e: escape speed
        f: galactic dm  halo vdf to transform
        """
       
        esc = esc_sun[0]
        sun = esc_sun[1]
        #n = integrate.quad(f, 0, v_e,args=(v_r,v_e))[0]
        result = integrate.quad(transform, -1, 1, args =(u, f, esc, sun))[0]
        return  u * result

def get_n(f):
    def intern(u):
        return gal_to_sun(u,f)
    n= 1/quad(intern, 0., v_esc + v_Sun)[0]
    return n


################# Alternative functions (ALT) for alternative purpuse
#def transformALT(x, u, f, es, sun, nor):
#        """
#        this is the galilean transformations that corresponds to the
#        expression 2.2  in the CHOI paper; the part inside the integral
#        """
#        v_S = v_Sun
#        if sun !=0:
#                v_S = sun
#        v = sqrt((u**2) + (v_S**2) + (2 * u * v_S * x)) 
#        return f(v,transform = True, esc=es, n = nor)
#
#def gal_to_sunALT(u, f, n, esc_sun=[v_esc, v_Sun]): 
#        """ 
#        This is the transformation of a generic VDF in the galactic halo frame of 
#        reference to the Sun f.o.r where
#        v_S: speed of the Sun
#        v_r: dispersion speed of the halo model
#        v_e: escape speed
#        f: galactic dm  halo vdf to transform
#        """
#       
#        esc = esc_sun[0]
#        sun = esc_sun[1]
#        #n = integrate.quad(f, 0, v_e,args=(v_r,v_e))[0]
#        result = integrate.quad(transformALT, -1, 1, args =(u, f, esc, sun,n))[0]
#        return  u * result




############################################################################
###############  in the galatic halo frame of reference  ###################

def choi_sun(u, v_S, v_r, v_e):
        """ 
        Choi et al maxwellian VDF 
        in the frame of reference of the sun
        """
        f = lambda v: (v**2) * (1/(v_r**3)) * exp((-3. * (v**2)) / (2. * v_r**2 ))
        n = integrate.quad(f, 0, v_e)[0]

        def transform(x, u, v_S, v_r, v_e):
   
                v = sqrt((u**2) + (v_S**2) + (2 * u * v_S * x)) 
                con =(4./sqrt(np.pi)) * (3./2.)**(3./2.)
                if v > v_e:
                        return 1e-30
                else:
                        t = (1/(v_r**3)) * exp((-3. * (v**2)) / (2. * v_r**2 ))
                return t
        result = integrate.quad(transform, -1, 1, args =(u, v_S, v_r, v_e))[0]
        return (1/n) * u * result

def mao_sun(u, v_S, v_r, v_e):
        """
        Mao et all empiric expression
        in frame of reference of the Sun
        """
        p = 7.2
        f = lambda v: (v**2) * exp(-v/v_r) * (((v_e**2) - (v**2))**p)
        n = integrate.quad(f, 0, v_e)[0]

        def transform(x, u, v_S, v_r, v_e):
                v = sqrt((u**2) + (v_S**2) + (2 * u * v_S * x))
                if v > v_e:
                        return 1e-60
                else:
                        t1 = exp(-v/v_r)
                        t2 = ((v_e**2) - (v**2))**p
                        return t1*t2
        result = integrate.quad(transform, -1, 1, args =(u, v_S, v_r, v_e))[0]
        return (1/n) * u * result

def lin_sun(u, v_S, v_r, v_e):
        """
        Tsallis expression in the
        Sun frame of reference
        from ling et al.
        """
        q = 0.77
        p = q / (1 - q)
        f = lambda v: (v**2) * ((1-((1-q) * (v**2 / v_r**2)))**p)
        n = integrate.quad(f, 0, v_e)[0]

        def transform(x, u, v_S, v_r, v_e):
                v = sqrt((u**2) + (v_S**2) + (2 * u * v_S * x))
                if v > v_e:
                        return 1e-30
                else:
                        t = (1-((1-q) * (v**2 / v_r**2)))
                        return t**p
        result = integrate.quad(transform, -1, 1, args =(u, v_S, v_r, v_e))[0]
        return (1/n) * u * result


####################### Bozorgnia's Paper ################################
#################### halo frame of reference #############################


def general_max(v, args):
    """
    Generalized maxwellian from Bozorgnia paper
    :param v: velocity modulus
    :param args: array of two required arguments (v_r, alpha)
           v_r: dispertion velocity
           alpha: power
    :return: value of the de distriburion depending on the v
    """
    global v_esc
    v_r = args[0]
    alpha = args[1]
    f = lambda v:(v ** 2) * exp(-(v / v_r) ** (2 * alpha))
    n = integrate.quad(f, 0, v_esc)[0]
    if v > v_esc:
        return 1e-40
    else:
        return (1 / n) * (v ** 2) * exp(-(v / v_r) ** (2 * alpha))


def maxwellian(v, v_r):
    """
    Standard maxwellian from Bozorgnia paper
    :param v: velocity modulus
    :param args: in this case alpha = 1
    so it gives only the dispession velocity
    :return: value of the de distriburion depending on the v
    """
    global v_esc
    f = lambda v: (v ** 2) * exp(-(v / v_r)**2)
    n = integrate.quad(f, 0, v_esc)[0]
    if v > v_esc:
        return 1e-40
    else:
        return (1 / n) * (v ** 2) * exp(-(v / v_r)**2)


def mao_boz_gal(v, args):
    """
    Mao expression that will be used  for the Borzorgnia
    simulations.
    :param v: velocity modulus
    :param args: (v_r,p)
    :return:
    """
    global v_esc
    v_r = args[0]
    p = args[1]
    f = lambda v: (v**2) * exp(-v / v_r) * (((v_esc ** 2) - (v ** 2)) ** p)
    n = integrate.quad(f, 0, v_esc)[0]
    if v > v_esc:
        return 1e-40
    else:
        t1 = exp(-v / v_r)
        t2 = ((v_esc ** 2) - (v ** 2)) ** p
        return (1/n) * v**2 * t1 * t2


def lisanti_gal(v,args):
    """
    lisanti distribution fit function from Borzorgnia
    :param v: velocity modulus
    :param args: (v_r, k )
    :return:
    """
    global v_esc
    v_r = args[0]
    k = args[1]
    f = lambda v: (v**2) * ((exp(((v_esc**2)-(v**2)) / (k * v_r**2))-1)**k)
    n = integrate.quad(f, 0, v_esc)[0]
    if v > v_esc:
        return 1e-40
    else:
        t1 = (((v_esc**2)-(v**2)) / (k * v_r**2))
        t2 = ((exp(t1)-1)**k)
        return (1/n) * v**2 * t2


#################### SUN frame of reference #############################


def general_max_Sun(u, args):
    """
    In the frame of reference of the Sun
    Generalized maxwellian from Bozorgnia paper
    :param v: velocity modulus
    :param args: array of two required arguments (v_r, alpha)
           v_r: dispertion velocity
           alpha: power
    :return: value of the de distriburion depending on the v
    """
    global v_esc
    global v_Sun
    v_r = args[0]
    alpha = args[1]

    def transform(x, u, v_Sun, v_r, v_esc):
        v = sqrt((u ** 2) + (v_Sun ** 2) + (2 * u * v_Sun * x))
        f = lambda v:(v ** 2) * exp(-(v / v_r) ** (2 * alpha))
        n = integrate.quad(f, 0, v_esc)[0]
        if v > v_esc:
            return 1e-40
        else:
            return (1/n) * exp(-(v / v_r) ** (2 * alpha))
    result = integrate.quad(transform, -1, 1, args=(u, v_Sun, v_r, v_esc))[0]
    return u * result

def maxwellian_Sun(u, v_r):
    """
    In the frame of reference of the Sun
    Standard maxwellian from Bozorgnia paper
    :param v: velocity modulus
    :param args: in this case alpha = 1
    so it gives only the dispession velocity
    :return: value of the de distriburion depending on the v
    """
    global v_esc
    global v_Sun

    def transform(x, u, v_Sun, v_r, v_esc):
        v = sqrt((u ** 2) + (v_Sun ** 2) + (2 * u * v_Sun * x))
        f = lambda v: (v ** 2) * exp(-(v / v_r)**2)
        n = integrate.quad(f, 0, v_esc)[0]
        if v > v_esc:
            return 1e-40
        else:
            return (1/n) * exp(-(v / v_r)**2)
    result = integrate.quad(transform, -1, 1, args=(u, v_Sun, v_r, v_esc))[0]
    return u * result


def mao_boz_gal_Sun(u, args):
    """
    Frame of reference of the Sun
    Mao expression that will be used  for the Borzorgnia
    simulations.
    :param v: velocity modulus
    :param args: (v_r,p)
    :return:
    """
    global v_esc
    global v_Sun
    v_r = args[0]
    p = args[1]

    def transform(x, u, v_Sun, v_r, v_esc):
        v = sqrt((u ** 2) + (v_Sun ** 2) + (2 * u * v_Sun * x))
        f = lambda v: (v**2) * exp(-v / v_r) * (((v_esc ** 2) - (v ** 2)) ** p)
        n = integrate.quad(f, 0, v_esc)[0]
        if v > v_esc:
            return 1e-40
        else:
            t1 = exp(-v / v_r)
            t2 = ((v_esc ** 2) - (v ** 2)) ** p
            return (1/n) * t1 * t2

    result = integrate.quad(transform, -1, 1, args=(u, v_Sun, v_r, v_esc))[0]
    return u * result


def lisanti_Sun(u,args):
    """
    Frame of reference of the Sun
    lisanti distribution fit function from Borzorgnia
    :param v: velocity modulus
    :param args: (v_r, k )
    :return:
    """
    global v_esc
    global v_Sun
    v_r = args[0]
    k = args[1]

    def transform(x, u, v_Sun, v_r, v_esc):
        v = sqrt((u ** 2) + (v_Sun ** 2) + (2 * u * v_Sun * x))
        f = lambda v: (v**2) * ((exp(((v_esc**2)-(v**2)) / (k * v_r**2))-1)**k)
        n = integrate.quad(f, 0, v_esc)[0]
        if v > v_esc:
            return 1e-40
        else:
            t1 = (((v_esc**2)-(v**2)) / (k * v_r**2))
            t2 = ((exp(t1)-1)**k)
            return  (1/n) * t2
    result = integrate.quad(transform, -1, 1, args=(u, v_Sun, v_r, v_esc))[0]
    return u * result
#############################################################


def get_disc_fdu(function, n_bins, histogram=False):
    """
    gets a discrete fdu in the sun frame of reference
    from de selected function
    and with a defined number of bins n_bins
    """
    bins = []
    hist = []
    prob = []
    tmp = []
    width = 700/n_bins
    for i in np.arange(0,700,1):
             if i % width == 0:
                     bins.append(i)
                     tmp.append(gal_to_sun(i, function))
    a = sum(tmp) 
    prob = [x / a  for x in tmp]
    if histogram == False:
            return bins, prob
    else:
            hist = [x / (a * width)  for x in tmp]
            return bins, prob, hist


######################################################
################FROM CLUSTER##########################
######################################################

def get_center_hist(x,y,z):
    binsize = 0.2
    cutoff= -12.8/binsize
    v_array = np.arange(0,np.max(x),binsize)
    histX, bins = np.histogram(x,bins=v_array)
    x_center = bins[np.where(histX==np.max(histX[:cutoff]))][0]
    v_array = np.arange(0,np.max(y),binsize)
    histY, bins = np.histogram(y,bins=v_array)
    y_center = bins[np.where(histY==np.max(histY[:cutoff]))][0]
    v_array = np.arange(0,np.max(z),binsize)
    histZ, bins = np.histogram(z,bins=v_array)
    z_center = bins[np.where(histZ==np.max(histZ[:cutoff]))][0]
    return x_center,y_center,z_center

def get_center_dens(dens_file):
    uns = CunsIn(dens_file,"all","all",False)
    ok = uns.nextFrame("")
    ok,rho = uns.getArrayF(a,"rho")
    ok,pos = uns.getArrayF(a,"pos")
    uns.close()
    rmi = np.where(rho==np.max(rho))[0][0]*3
    return pos[rmi],pos[rmi+1],pos[rmi+2]
    

def get_fdv(file,lowbound,upbound,verbose=False, unittokpc = 1):
    limmax,limmin = float(upbound),float(lowbound)
    a="all"
    uns = CunsIn(file,"halo",a,False)
    bits=""         # select properties, "" means all
    ok=uns.nextFrame(bits)   # load data from disk
    ok,vel = uns.getArrayF(a,"vel")
    ok,pos = uns.getArrayF(a,"pos")
    uns.close()
    vel = vel*unittokpc
    x = pos[0::3]
    y = pos[1::3]
    z = pos[2::3]
    center = get_center(x,y,z)
    x = x-center[0]
    y = y-center[1]
    z = z-center[2]
    r = np.sqrt((x**2)+(y**2)+(z**2))
    
    vel_modulus = np.sqrt((vel[0::3]**2)+(vel[1::3]**2)+(vel[2::3]**2))
    r_tmp = r[np.where(r<limmax)]
    v_modtmp = vel_modulus[np.where(r < limmax )]
    
    v_mod = v_modtmp[np.where(r_tmp>limmin)]
    bin_size = .1*unittokpc
    v_array = np.arange(0,np.max(v_mod),bin_size)
    hist, bins = np.histogram(v_mod,bins=v_array,normed=1)
    if verbose==True:
        print file[-15:]
        print ok
        print "center =", center
        print "n_part with "+str(limmin)+" < r < "+str(limmax)+" =",len(v_mod),
        print "total number of particles =",len(vel_modulus)
    
    return bins[:-1],hist,file[-15:]



def get_fdv_disc(files,lowboundr,upboundr,thickness,**kwargs):
    binsize = kwargs.get('binsize',1)
    comp = kwargs.get('comp',"halo")
    dens_file = kwargs.get('dens_file',None)
    density = kwargs.get('density',False)
    disk = kwargs.get('disk',False)
    plane = kwargs.get('plane',"z")
    theta = kwargs.get('theta',45.0)
    rotate = kwargs.get('rotate',False)
    verbose = kwargs.get('verbose',False)
    simutoKPC =kwargs.get('simutoKPC',1)
    simutoKMS =kwargs.get('simutoKMS',1)
    ##################################### 
    limmax,limmin = float(upboundr),float(lowboundr)
    a="all"
    uns = CunsIn(files,comp,a,False)
    bits=""         # select properties, "" means all
    ok=uns.nextFrame(bits)   # load data from disk
    ok,vel = uns.getArrayF(a,"vel")
    ok,pos = uns.getArrayF(a,"pos")
    uns.close()
    ####################################
    vel = vel*simutoKMS
    pos = pos*simutoKPC
    x = pos[0::3]
    y = pos[1::3]
    z = pos[2::3]
    if density == False:
        center = get_center_hist(x,y,z)
    else:
        center = get_center_dens(dens_file)
    ######### recenter box #############    
    x = x-center[0]
    y = y- center[1]
    z = z-center[2]
    ######## rotating part #############
    if rotate==True:
        theta = np.radians(theta)
        x = x / np.cos(theta)
        z = z / np.sin(theta)
    ########   r and |v|   #############
    r = np.sqrt((x**2)+(y**2)+(z**2))
    vel_modulus = np.sqrt((vel[0::3]**2)+(vel[1::3]**2)+(vel[2::3]**2))
    ####################################
    ######  shell particles  ###########
    if disk == False:
        r = r[(r<=limmax)&(r>=limmin)]
        v = vel_modulus[(r<=limmax)&(r>=limmin)]
    ######  disc particles  ###########
    else:
        if plane=="z":
            condition  = (r<=limmax)&(r>=limmin)&(abs(z)<=thickness)
            r = r[condition]
            v = vel_modulus[condition]
        if plane=="x": 
            condition  = (r<=limmax)&(r>=limmin)&(abs(x)<=thickness)
            r = r[condition]
            v = vel_modulus[condition]
        if plane=="y": 
            condition  = (r<=limmax)&(r>=limmin)&(abs(y)<=thickness)
            r = r[condition]
            v = vel_modulus[condition]
    ####### histogram data ############
    bin_size = binsize
    v_array = np.arange(0,np.max(v),bin_size)
    hist, bins = np.histogram(v,bins=v_array, normed=1)
    ######  verbose printing ##########
    if verbose==True:
        print "###########################"
        print files[-15:]
        print ok
        print "total halo particles =", len(pos)/3
        print "center =", center
        if disk==False:
            print "n_part in shell with "+str(limmin)+" < r < "+str(limmax)+" =",len(v)
        else:
            print "n_part in disc with "+str(limmin)+" < r < "+str(limmax)+" and |"+plane+"| < "+str(thickness)+"=",len(v)

    return bins[:-1],hist,files[-15:]


class info:
    def __init__(self,path):
        vars=dict()
        with open(glob.glob(path+'/info_?????'+'.txt')[0]) as outinfo:
            for line in outinfo:
                eq_index = line.find('=')
                if eq_index == -1:
                    continue
                var_name = line[:eq_index].strip()
                number = float(line[eq_index + 1:].strip())
                vars[var_name] = number
                if var_name == "unit_t":
                    break   
        aexp=vars["aexp"]
	self.z = (1-aexp)/aexp
        self.time = vars["time"]*aexp
        self.Msunkg = 1.99844*10**30
        self.pctocm = 3.08567758*10**18
        self.cmtopc = 1/self.pctocm
        self.unitl=vars["unit_l"]/(3.08*10**18)*self.pctocm
        self.unitd=vars["unit_d"]/(self.pctocm/(3.08*10**18))**3
        self.unitt=vars["unit_t"]
        self.simutokpc=self.unitl/self.pctocm/1000
        self.simutoMsun=(self.unitd*self.unitl**3)/1000/self.Msunkg
        self.unitsimutoMsunkpc3=self.unitd*self.pctocm**3/1000/self.Msunkg
        self.unitsimutokms = self.unitl/10**5/self.unitt

    


class Halo:
    def __init__(self, file_path, component,**kwargs):
        simutokpc = kwargs.get('simutokpc',1.)
        simutokms = kwargs.get('simutokms',1.)
        hsml = kwargs.get('hsml',False)
        simutoMsun = kwargs.get('simutoMsun', 1.)
        dens = kwargs.get('dens',False)        
        center = kwargs.get('center',[0,0,0])
        halo_vel = kwargs.get('halo_vel',[0.,0.,0.])
        self.comp = component
        msuntokg = 1.989e30 
        kgtoGeV = 1/1.783e-27
        kpctocm = 3.086e21
        simutoGeVcm3 = (simutoMsun*msuntokg*kgtoGeV) / (simutokpc*kpctocm)**3 
        uns = CunsIn(file_path,component,"all",False)
        if uns.isValid()!=True:
            sys.exit("\n\n\n\n\n\n\nERROR:"+file_path+" is not a valid file !!!!!\n\n\n\n\n")
        ok=uns.nextFrame("") 
        ok,pos = uns.getArrayF("all","pos")
        ok, vel = uns.getArrayF("all","vel")
        if dens ==True:
            ok, rho = uns.getArrayF("all","rho")
            self.rho =  rho * simutoGeVcm3
        if hsml ==True:
            ok, self.hsml = uns.getArrayF("all","hsml")
        uns.close()
        
        ### coordinates ###
        pos = pos * simutokpc
        self.x,self.y,self.z = pos[0::3]-center[0],pos[1::3]-center[1],pos[2::3]-center[2]
        self.R = np.sqrt((self.x**2)+(self.y**2))
        self.r = np.sqrt((self.x**2)+(self.y**2)+(self.z**2))
        ### velocities ###
        vel = vel * simutokms
        self.vx,self.vy,self.vz = vel[0::3] - halo_vel[0],vel[1::3] - halo_vel[1],vel[2::3] - halo_vel[2]
        self.v = np.sqrt((self.vx**2) + (self.vy**2) + (self.vz**2))
        self.vR = (self.vx*self.x + self.vy*self.y)/ self.R
        self.vr = (self.vx*self.x + self.vy*self.y + self.vz*self.z)/ self.r
        self.vphi = (-self.vx*self.y + self.vy*self.x )/ self.R
        ### densities ### 
        
    def get_shell(self,lim_min,lim_max):
        shell_cond = (self.r<=lim_max)&(self.r>=lim_min)
        self.v_shell = self.v[shell_cond]
        self.vr_shell = self.vr[shell_cond]
        self.vR_shell = self.vR[shell_cond]
        self.vphi_shell = self.vphi[shell_cond]
        self.vz_shell = self.vz[shell_cond]
    
    def get_ring_xy(self,lim_min,lim_max, thickness):
        ringXY_cond = (self.r<=lim_max)&(self.r>=lim_min)&(abs(self.z)<=thickness)
        self.v_ringXY = self.v[ringXY_cond]
        self.vr_ringXY = self.vr[ringXY_cond]
        self.vR_ringXY = self.vR[ringXY_cond]
        self.vphi_ringXY = self.vphi[ringXY_cond]
        self.vz_ringXY = self.vz[ringXY_cond]
        self.x_ringXY = self.x[ringXY_cond]
        self.y_ringXY = self.y[ringXY_cond]
        self.z_ringXY = self.z[ringXY_cond]
        self.vx_ringXY = self.vx[ringXY_cond]
        self.vy_ringXY = self.vy[ringXY_cond]
        
    def get_ring_yz(self,lim_min,lim_max, thickness):
        ringYZ_cond = (self.r<=lim_max)&(self.r>=lim_min)&(abs(self.x)<=thickness)
        self.v_ringYZ = self.v[ringYZ_cond]
        self.vr_ringYZ = self.vr[ringYZ_cond]
        self.vR_ringYZ = self.vR[ringYZ_cond]
        self.vphi_ringYZ = self.vphi[ringYZ_cond]
        self.vz_ringYZ = self.vz[ringYZ_cond]
        
    def get_ring_zx(self,lim_min,lim_max, thickness):
        ringZX_cond = (self.r<=lim_max)&(self.r>=lim_min)&(abs(self.y)<=thickness)
        self.v_ringZX = self.v[ringZX_cond]
        self.vr_ringZX = self.vr[ringZX_cond]
        self.vR_ringZX = self.vR[ringZX_cond]
        self.vphi_ringZX = self.vphi[ringZX_cond]
        self.vz_ringZX = self.vz[ringZX_cond]
    
    def get_dens(self,lim_min,lim_max, thickness,plane="x"):
        shell_cond = (self.r<=lim_max)&(self.r>=lim_min)
        if plane == "x":
            plane_cond = (abs(self.z)<=thickness)
        elif plane == "y":
            plane_cond = (abs(self.x)<=thickness)
        elif plane == "z":
            plane_cond = (abs(self.y)<=thickness)
        else:
            plane_cond = (True)
        condition = shell_cond & plane_cond
        self.rho_ring = self.rho[condition]
        self.rho_shell = self.rho[shell_cond]

  

############################



class Fit:
    def __init__(self,ndim,**var):
        self.range1 = var.get('range1',[0.,10.])
        self.range2 = var.get('range2',[0.,10.])
        self.range3 = var.get('range3',[0.,10.])
        self.range4 = var.get('range4',[0.,10.])
        self.ndim = ndim
        if self.ndim == 1:
            self.pos_min = np.array([100.0])
            self.pos_max = np.array([300.])
            self.condition = ()
        elif self.ndim == 2: 
            self.pos_min = np.array([-5.0,0.])
            self.pos_max = np.array([5.0,10.])
        elif self.ndim == 3: 
            self.pos_min = np.array([-20.0,0.,-2.])
            self.pos_max = np.array([20.0,600.,2.])
        elif self.ndim == 4: 
            self.pos_min = np.array([.0,-300.,-2.,-300.])
            self.pos_max = np.array([1.,300.,2.,300])
        else:
            print "ERROR: ndim > 4 not supported"
            sys.exit()
            
    def bestfit(self,func,x,y):
        yerr = 0.002 + 0.005*np.random.rand(len(x))
        # Now, let's setup some parameters that define the MCMC
        self.nwalkers = 500

        # Initialize the chain
        # Choice 1: chain uniformly distributed in the range of the parameters
        
        self.psize = self.pos_max - self.pos_min
        self.pos = [self.pos_min + self.psize*np.random.rand(self.ndim) for i in range(self.nwalkers)]

        # As prior, we assume an 'uniform' prior (i.e. constant prob. density)
        if self.ndim==1:
            def lnprior(theta):
                v0 = theta
                if self.range1[0] < v0 < self.range1[1]: 
                    return 0.0
                return -np.inf
            
            def lnlike(theta, x, y, yerr):
                v0 = theta
                alpha = 1.
                model = func(x,v0)
                return -0.05*(np.sum( ((y-model)/yerr)**2. ))
            
        if self.ndim==2:
            def lnprior(theta):
                mu,v0 = theta
                if self.range1[0] < mu < self.range1[1] and  \
                   self.range2[0] < v0 < self.range2[1]:
                    return 0.0
                return -np.inf
            
            def lnlike(theta, x, y, yerr):
                mu, v0 = theta
                alpha = 1.
                model = func(x,mu,v0)
                return -0.05*(np.sum( ((y-model)/yerr)**2. ))
            
        if self.ndim==3:
            def lnprior(theta):
                mu,v0,alpha = theta
                if self.range1[0] < mu < self.range1[1] and  \
                    self.range2[0] < v0 < self.range2[1] and  \
                     self.range3[0] < alpha < self.range3[1]:
                            return 0.0
                return -np.inf
            
            def lnlike(theta, x, y, yerr):
                mu, v0, alpha = theta
                alpha = 1.
                model = func(x,mu,v0, alpha)
                return -0.05*(np.sum( ((y-model)/yerr)**2. ))
            
        if self.ndim==4:
            def lnprior(theta):
                f,mu,v0,alpha = theta
                if self.range1[0] < f < self.range1[1] and  \
                    self.range2[0] < mu < self.range2[1] and  \
                     self.range3[0] < v0 < self.range3[1] and \
                        self.range4[0] < alpha < self.range4[1]:
                            return 0.0
                return -np.inf
           
            def lnlike(theta, x, y, yerr):
                frac,v01, mu, v02 = theta
                alpha = 1.
                model = func(x,frac,v01,mu,v02)
                return -0.05*(np.sum( ((y-model)/yerr)**2. ))
                
        # As likelihood, we assume the chi-square. Note: we do not even need to normalize it.
        

        def lnprob(theta, x, y, yerr):
            lp = lnprior(theta)
            if not np.isfinite(lp):
                return -np.inf
            return lp + lnlike(theta, x, y, yerr)
        print self.ndim
        sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim,
                                lnprob, args=(x, y, yerr))
        self.pos, self.prob, self.state  = sampler.run_mcmc(self.pos, 300)
        sampler.reset()
        self.pos, self.prob, self.state  = sampler.run_mcmc(self.pos, 1000)
        samples = sampler.flatchain
        samples.shape
        print len(samples)
        sam=samples[-100:]
        if self.ndim==1:
            self.expected = np.array([func(x,i[0]) for i in sam ])
        if self.ndim==2:
            self.expected = np.array([func(x,i[0],i[1]) for i in sam ])
        if self.ndim==3:
            self.expected = np.array([func(x,i[0],i[1],i[2]) for i in sam ])
        if self.ndim==4:
            self.expected = np.array([func(x,i[0],i[1],i[2],i[3]) for i in sam ])
        self.observed = np.array([y for i in sam])
        chi_tmp = ((self.observed-self.expected)**2)/self.expected
        chi2 = [np.sum(i) for i in chi_tmp]
        print chi2#np.is
        print np.nanmin(chi2), np.where(chi2==np.nanmin(chi2))[0][0]
        self.params = sam[np.where(chi2==np.nanmin(chi2))[0][0]]
        return self.params, np.nanmin(chi2)
    
    
   
