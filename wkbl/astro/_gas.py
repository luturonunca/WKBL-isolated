"""
gas  and AMR routine
"""
import sys
import math
import glob
import cmath
import subprocess
import numpy as np
from py_unsio import *
import scipy.special as sp
from numpy import exp, sqrt
import nbody_essentials as nbe
import scipy.integrate as integrate
from sklearn.neighbors import KDTree

class _gas:
    def __init__(self, file_path,p, **kwargs):
        self._p = p
        self._center_history = np.array([0.,0.,0.]) 
        self.uns = CunsIn(file_path,"gas","all",False)
        dens = kwargs.get('dens',False)
        self.get_sigma = kwargs.get('virial',False)
        comov = kwargs.get('comov',False)
        self.halo_vel = kwargs.get('halo_vel',[0.,0.,0.])    ##########

        if self.uns.isValid()!=True:
            sys.exit("\n\n\n\n\n\n\nERROR:"+file_path+" is not a valid file !!!!!\n\n\n\n\n")
        ok = self.uns.nextFrame("")
        ok, pos = self.uns.getArrayF("gas","pos")
        ok, vel = self.uns.getArrayF("gas","vel")
        ok, mass = self.uns.getArrayF("gas","mass")
        ok, self.metal = self.uns.getArrayF("gas","metal")
        ok, temp = self.uns.getArrayF("gas","temp")
        ok, temp2 = self.uns.getArrayF("hydro","4")
        ok, pot = self.uns.getArrayF("gas","pot")
        ok, self.id = self.uns.getArrayI("gas","id")
        ok, rho = self.uns.getArrayF("gas","rho")
        ok, hsml =  self.uns.getArrayF("gas","hsml")
        if (self.get_sigma):
            ok, sigma = self.uns.getArrayF("hydro","7")
            self.sigma2 = sigma*self._p.simutokpc 
        
        ### coordinates ###
        vel = vel * self._p.simutokms
        if (comov):
            pos = pos * self._p.simutokpc / self._p.aexp
        else:
            pos = pos * self._p.simutokpc
        gamma = 1.666
        self.temp = (gamma-1.0) * temp * self._p.simutoKelvin #Kelvin
        self.rho =  rho * self._p.simutoMsun * (self._p.simutokpc)**-3
        self.temp2 = temp2 * self._p.simutoKelvin #Kelvin
        self.pot = pot #* self._p.simutokms**2 
        self.pos3d = pos.reshape(len(pos)/3,3)
        self.center_rho_max = self.pos3d[np.where(self.rho==self.rho.max())][0]
        self.vel3d = vel.reshape(len(vel)/3,3)
        self.mass = mass * self._p.simutoMsun
        self.hsml = hsml * self._p.simutokpc

    def halo_Only(self, center, n, r200):
        in_halo = nbe.all_inside(self.pos3d, center, n*r200)
        self.pos3d = self.pos3d[in_halo] - center
        self.mass = self.mass[in_halo]
        if (self.get_sigma):
            self.sigma2 = self.sigma2[in_halo]
        self.hsml = self.hsml[in_halo]
        self.metal = self.metal[in_halo]
        self.temp = self.temp[in_halo]
        self.temp2 = self.temp2[in_halo]
        self.pot = self.pot[in_halo]
        self.rho = self.rho[in_halo]
        self.vel3d = self.vel3d[in_halo]
        self.id = self.id[in_halo]
        self.R = np.sqrt((self.pos3d[:,0]**2)+(self.pos3d[:,1]**2))
        self.r = np.sqrt((self.pos3d[:,0]**2)+(self.pos3d[:,1]**2)+(self.pos3d[:,2]**2))
        ### velocities ###
        average_v = np.array([np.mean(self.vel3d[:,0]),np.mean(self.vel3d[:,1]),np.mean(self.vel3d[:,2])])
        self.vel3d = self.vel3d - average_v
        vx,vy,vz = self.vel3d[:,0],self.vel3d[:,1],self.vel3d[:,2]
        self.v = np.sqrt((vx**2) + (vy**2) + (vz**2))
        self.vR = (vx*self.pos3d[:,0] + vy*self.pos3d[:,1])/ self.R
        self.vr = (vx*self.pos3d[:,0] + vy*self.pos3d[:,1] + vz*self.pos3d[:,2])/ self.r
        self.vphi = (-vx*self.pos3d[:,1] + vy*self.pos3d[:,0] )/ self.R
        self.vtheta = (self.vR*self.pos3d[:,2] - vz*self.R) / self.r
        #### other params ###
        self.total_m =  np.sum(self.mass)

    def rotate(self,T):
        pos = self.pos3d
        self.pos3d = nbe.matrix_vs_vector(T,pos)

    def shift(self,center):
        self.pos3d = self.pos3d - center
        self._center_history = np.vstack((self._center_history,center))

