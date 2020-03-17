#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

fit class to fit velocity distribution funcitons
on n-body data from ramses
f(u)/u
for dark halos
author : arturonunez25@gmail.com

---------- FITTIN VDF SECTION -------------
"""
import sys
import math
import glob
import cmath
<<<<<<< HEAD
#import emcee
=======
>>>>>>> python3Ver
import subprocess
import numpy as np
from unsio import *
import scipy.special as sp
from numpy import exp, sqrt
import scipy.integrate as integrate
from wkbl.inputing import dm_input
quad = integrate.quad
sys.dont_write_bytecode = True


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
        self.params = sam[np.where(chi2==np.nanmin(chi2))[0][0]]
        return self.params, np.nanmin(chi2)
