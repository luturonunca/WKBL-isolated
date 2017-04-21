#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

fit class to fit velocity distribution funcitons
on nbody data from ramses
f(u)/u
for dark halos
author : arturonunez25@gmail.com

---------- FITTIN VDF SECTION -------------
"""
import sys
import math
import glob
import cmath
import emcee
import subprocess
import numpy as np
from py_unsio import *
import scipy.special as sp
from numpy import exp, sqrt
from scipy.special import gamma
from wkbl.inputing import dm_input
import scipy.integrate as integrate
quad = integrate.quad
sys.dont_write_bytecode = True

class DG_fit:
    """
    performs a fit of a double gaussian
    """
    def __init__(self,**var):
        self.pos_min = var.get('mins',np.array([0,0.,-400.,0.]))
        self.pos_max = var.get('maxs',np.array([.5,600.,400.,600]))

    def bestfit(self,x,y):
        def myGaussian(x,mu,v0):
            alpha = 1.
            n = 2 * v0 * gamma(1+(1/(2*alpha)))
            return (1/n) * np.exp(-((x-mu)**2/(v0**2))**alpha)
        
        def myDoubleGaussianFunc(x,frac,v01,mu2,v02):
            return frac*myGaussian(x,mu2,v01)+(1-frac)*myGaussian(x,0,v02)
        self.ndim = 4
        func = myDoubleGaussianFunc
        yerr = 0.000002 + 0.00005*np.random.rand(len(x))
        # Now, let's setup some parameters that define the MCMC
        self.nwalkers = 500
        # Initialize the chain
        # Choice 1: chain uniformly distributed in the range of the parameters
        self.psize = self.pos_max - self.pos_min
        self.pos = [self.pos_min + self.psize*np.random.rand(self.ndim) for i in range(self.nwalkers)]
        # As prior, we assume an 'uniform' prior (i.e. constant prob. density)
        def lnprior(theta):
            f,mu,v0,alpha = theta
            if self.pos_min[0] < f < self.pos_max[0] and  \
                   self.pos_min[1] < mu < self.pos_max[1] and  \
                     self.pos_min[2] < v0 < self.pos_max[2] and \
                        self.pos_min[3] < alpha < self.pos_max[3]:
                            return 0.0
            return -np.inf

        def lnlike(theta, x, y,yerr):
            frac,v01, mu, v02 = theta
            alpha = 1.
            model = func(x,frac,v01,mu,v02)
            return -0.005*(np.sum( ((y-model)/yerr)**2. ))

        # As likelihood, we assume the chi-square. Note: we do not even need to normalize it.

        def lnprob(theta, x, y, yerr):
            lp = lnprior(theta)
            if not np.isfinite(lp):
                return -np.inf
            return lp + lnlike(theta, x, y, yerr)        
        sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim,
                                lnprob, args=(x, y, yerr))
        sampler.reset()
        self.pos, self.prob, self.state  = sampler.run_mcmc(self.pos, 300)
        sampler.reset()
        self.pos, self.prob, self.state  = sampler.run_mcmc(self.pos, 700)
        samples = sampler.flatchain
        samples.shape
        sam=samples[-100:]
        self.expected = np.array([func(x,i[0],i[1],i[2],i[3]) for i in sam ])
        self.observed = np.array([y for i in sam])
        chi_tmp = ((self.observed-self.expected)**2)/self.expected
        chi2 = np.array([np.sum(i) for i in chi_tmp])
        chi2 = chi2[np.where(chi2>0)]
        self.params = sam[np.where(chi2==np.nanmin(chi2))[0][0]]
        return self.params, np.nanmin(chi2)


class SG_fit:
    """
    performs a fit of a Single gaussian
    """
    def __init__(self,**var):
        self.pos_min = var.get('mins',np.array([-500,0.]))
        self.pos_max = var.get('maxs',np.array([500,600.]))

    def bestfit(self,x,y):
        def myGaussian(x,mu,v0):
            alpha = 1.
            n = 2 * v0 * gamma(1+(1/(2*alpha)))
            return (1/n) * np.exp(-((x-mu)**2/(v0**2))**alpha)
        
        self.ndim = 2
        func = myGaussian
        yerr = 0.000002 + 0.00005*np.random.rand(len(x))
        # Now, let's setup some parameters that define the MCMC
        self.nwalkers = 500
        # Initialize the chain
        # Choice 1: chain uniformly distributed in the range of the parameters
        self.psize = self.pos_max - self.pos_min
        self.pos = [self.pos_min + self.psize*np.random.rand(self.ndim) for i in range(self.nwalkers)]
        # As prior, we assume an 'uniform' prior (i.e. constant prob. density)
        def lnprior(theta):
            mu,v0= theta
            if self.pos_min[0] < mu < self.pos_max[0] and  \
                   self.pos_min[1] < v0 < self.pos_max[1]:
                            return 0.0
            return -np.inf

        def lnlike(theta, x, y,yerr):
            mu, v01 = theta
            alpha = 1.
            model = func(x,mu,v01)
            return -0.005*(np.sum( ((y-model)/yerr)**2. ))

        # As likelihood, we assume the chi-square. Note: we do not even need to normalize it.

        def lnprob(theta, x, y, yerr):
            lp = lnprior(theta)
            if not np.isfinite(lp):
                return -np.inf
            return lp + lnlike(theta, x, y, yerr)        
        sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim,
                                lnprob, args=(x, y, yerr))
        sampler.reset()
        self.pos, self.prob, self.state  = sampler.run_mcmc(self.pos, 300)
        sampler.reset()
        self.pos, self.prob, self.state  = sampler.run_mcmc(self.pos, 700)
        samples = sampler.flatchain
        samples.shape
        sam=samples[-100:]
        self.expected = np.array([func(x,i[0],i[1]) for i in sam ])
        self.observed = np.array([y for i in sam])
        chi_tmp = ((self.observed-self.expected)**2)/self.expected
        chi2 = np.array([np.sum(i) for i in chi_tmp])
        chi2 = chi2[np.where(chi2>0)]
        self.params = sam[np.where(chi2==np.nanmin(chi2))[0][0]]
        return self.params, np.nanmin(chi2)


class GM_fit:
    """
    performs a fit of a Generalize Maxwellian
    """
    def __init__(self,**var):
        self.pos_min = var.get('mins',np.array([0,-10.]))
        self.pos_max = var.get('maxs',np.array([500,10.]))

    def bestfit(self,x,y):
        def myGenMaxwellian(x,v0,alpha):
            n = v0**3 * gamma(1. + (3./(2.*alpha))) / 3.
            return (x**2/n) * np.exp(-(x/v0)**(2*alpha))
        
        self.ndim = 2
        func = myGenMaxwellian
        yerr = 0.000002 + 0.00005*np.random.rand(len(x))
        # Now, let's setup some parameters that define the MCMC
        self.nwalkers = 500
        # Initialize the chain
        # Choice 1: chain uniformly distributed in the range of the parameters
        self.psize = self.pos_max - self.pos_min
        self.pos = [self.pos_min + self.psize*np.random.rand(self.ndim) for i in range(self.nwalkers)]
        # As prior, we assume an 'uniform' prior (i.e. constant prob. density)
        def lnprior(theta):
            mu,v0= theta
            if self.pos_min[0] < mu < self.pos_max[0] and  \
                   self.pos_min[1] < v0 < self.pos_max[1]:
                            return 0.0
            return -np.inf

        def lnlike(theta, x, y,yerr):
            mu, v01 = theta
            alpha = 1.
            model = func(x,mu,v01)
            return -0.005*(np.sum( ((y-model)/yerr)**2. ))

        # As likelihood, we assume the chi-square. Note: we do not even need to normalize it.

        def lnprob(theta, x, y, yerr):
            lp = lnprior(theta)
            if not np.isfinite(lp):
                return -np.inf
            return lp + lnlike(theta, x, y, yerr)        
        sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim,
                                lnprob, args=(x, y, yerr))
        sampler.reset()
        self.pos, self.prob, self.state  = sampler.run_mcmc(self.pos, 300)
        sampler.reset()
        self.pos, self.prob, self.state  = sampler.run_mcmc(self.pos, 700)
        samples = sampler.flatchain
        samples.shape
        sam=samples[-100:]
        self.expected = np.array([func(x,i[0],i[1]) for i in sam ])
        self.observed = np.array([y for i in sam])
        chi_tmp = ((self.observed-self.expected)**2)/self.expected
        chi2 = np.array([np.sum(i) for i in chi_tmp])
        chi2 = chi2[np.where(chi2>0)]
        self.params = sam[np.where(chi2==np.nanmin(chi2))[0][0]]
        return self.params, np.nanmin(chi2)


class al_be_ga_fit:
    """
    performs a fit of a alpha beta gamma density profile fit
    """
    def __init__(self,**var):
        self.pos_min = var.get('mins',np.array([2e8,0.,-1.,-1.,-1.]))
        self.pos_max = var.get('maxs',np.array([2e10,10.,5.,5.,5.]))

    def bestfit(self,x,y,n):
        def abg_profile(x,p_s,r_s,al,be,ga):
            power =  (be - ga) / al
            denominator = ((x/r_s)**ga) * (1 + (x / r_s)**al)**power
            return p_s / denominator
        
        self.ndim = 5
        func = abg_profile
        yerr = 1 / np.sqrt(n)
        # Now, let's setup some parameters that define the MCMC
        self.nwalkers = 500
        # Initialize the chain
        # Choice 1: chain uniformly distributed in the range of the parameters
        self.psize = self.pos_max - self.pos_min
        self.pos = [self.pos_min + self.psize*np.random.rand(self.ndim) for i in range(self.nwalkers)]
        # As prior, we assume an 'uniform' prior (i.e. constant prob. density)
        def lnprior(theta):
            p_s, r_s, al, be, ga = theta
            if self.pos_min[0] < p_s < self.pos_max[0] and  \
                   self.pos_min[1] < r_s < self.pos_max[1] and  \
                     self.pos_min[2] < al < self.pos_max[2] and \
                        self.pos_min[3] < be < self.pos_max[3] and \
                          self.pos_min[4] < ga < self.pos_max[4]: 
                            return 0.0
            return -np.inf

        def lnlike(theta, x, y,yerr):
            p_s, r_s, al, be, ga = theta
            alpha = 1.
            model = func(x,p_s, r_s, al, be, ga)
            return -0.5*(np.sum( ((y-model)/yerr)**2. ))

        # As likelihood, we assume the chi-square. Note: we do not even need to normalize it.

        def lnprob(theta, x, y, yerr):
            lp = lnprior(theta)
            if not np.isfinite(lp):
                return -np.inf
            return lp + lnlike(theta, x, y, yerr) 
        
        sampler = emcee.EnsembleSampler(self.nwalkers, self.ndim,
                                lnprob, args=(x, y, yerr))
        sampler.reset()
        self.pos, self.prob, self.state  = sampler.run_mcmc(self.pos, 300)
        sampler.reset()
        self.pos, self.prob, self.state  = sampler.run_mcmc(self.pos, 700)
        samples = sampler.flatchain
        samples.shape
        sam=samples[-400:]
        self.expected = np.array([func(x,i[0],i[1],i[2],i[3],i[4]) for i in sam ])
        self.observed = np.array([y for i in sam])
        chi_tmp = ((np.log10(self.observed)-np.log10(self.expected))**2)/np.log10(self.expected)
        chi2 = np.array([np.sum(i) for i in chi_tmp])
        chi2 = chi2[np.where(chi2>0)]
        self.params = sam[np.where(chi2==np.nanmin(chi2))[0][0]]
        return self.params, np.nanmin(chi2)
