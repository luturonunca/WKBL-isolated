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
    def __init__(self, file_path,p,dens=True, **kwargs):
        self._p = p
        self._center_history = np.array([0.,0.,0.]) 
        self.uns = CunsIn(file_path,"gas","all",False)
        hsml = kwargs.get('hsml',False)
        dens = kwargs.get('dens',False)
        comov = kwargs.get('comov',False)
        self.halo_vel = kwargs.get('halo_vel',[0.,0.,0.])    ##########

        if self.uns.isValid()!=True:
            sys.exit("\n\n\n\n\n\n\nERROR:"+file_path+" is not a valid file !!!!!\n\n\n\n\n")
        ok = self.uns.nextFrame("")
        ok, pos = self.uns.getArrayF("all","pos")
        ok, vel = self.uns.getArrayF("all","vel")
        ok, mass = self.uns.getArrayF("all","mass")
        ok, temp = self.uns.getArrayF("gas","temp")
        ok, self.id = self.uns.getArrayI("all","id")
        ok, rho = self.uns.getArrayF("gas","rho")
        self.rho =  rho * self._p.simutoMsun / (self._p.simutokpc**3)
        ok, hsml = self.uns.getArrayF("gas","hsml")
        ### coordinates ###
        pos = pos * self._p.simutokpc
        vel = vel * self._p.simutokms
        if (comov):
            self.pos3d = pos.reshape(len(pos)/3,3) / self._p.aexp
        else:
            self.pos3d = pos.reshape(len(pos)/3,3)
        self.vel3d = vel.reshape(len(vel)/3,3)
        self.mass = mass * self._p.simutoMsun
        self.hsml = hsml * self._p.simutokpc
        self.center_rho_max = self.pos3d[np.where(self.rho == self.rho.max())]
        tokelvin =  self._p.mu * (self._p.simutokms**2) * self._p.cmtopc * 1e4 / self._p.kB # the 1e4 is to have km/s into kpc/s
        self.temp = temp * tokelvin 

    def halo_Only(self, center, n, r200):
        self.r = np.sqrt((self.pos3d[:,0]**2)+(self.pos3d[:,1]**2)+(self.pos3d[:,2]**2))
        in_halo = np.where(self.r <= n*r200)
        self.pos3d = self.pos3d[in_halo]
        self.mass = self.mass[in_halo]
        self.temp = self.temp[in_halo]
        self.hsml = self.hsml[in_halo]
        self.vel3d = self.vel3d[in_halo]
        self.id = self.id[in_halo]
        self.rho = self.rho[in_halo]
        self.R = np.sqrt((self.pos3d[:,0]**2)+(self.pos3d[:,1]**2))
        self.r = np.sqrt((self.pos3d[:,0]**2)+(self.pos3d[:,1]**2)+(self.pos3d[:,2]**2))
        self.phi = np.arctan2(np.copy(self.pos3d[:,1]),np.copy(self.pos3d[:,0]))
        ### velocities ###
        vx,vy,vz = self.vel3d[:,0],self.vel3d[:,1],self.vel3d[:,0]*self.pos3d[:,2]
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
        self.vel3d = nbe.matrix_vs_vector(T,self.vel3d)


    def shift(self,center):
        self.pos3d = self.pos3d - center
        self._center_history = np.vstack((self._center_history,center))

