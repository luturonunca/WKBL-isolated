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
        ok, pot = self.uns.getArrayF("gas","pot")
        ok, vel = self.uns.getArrayF("all","vel")
        ok, mass = self.uns.getArrayF("all","mass")
        ok, self.id = self.uns.getArrayI("all","id")
        ok, rho = self.uns.getArrayF("gas","rho")
        self.rho =  rho * self._p.simutoMsun / (self._p.simutokpc**3)
        ok, hsml = self.uns.getArrayF("gas","hsml")
        ### coordinates ###
        vel = vel * self._p.simutokms
        if (comov):
            pos = pos * self._p.simutokpc / self._p.aexp
        else:
            pos = pos * self._p.simutokpc
        
        self.pos3d = pos.reshape(len(pos)/3,3)
        self.vel3d = vel.reshape(len(vel)/3,3)
        self.mass = mass * self._p.simutoMsun
        self.hsml = hsml * self._p.simutokpc
        self.pot = pot * self._p.simutokms
        #self.center_com = nbe.real_center(self.pos3d,self.mass)

    def halo_Only(self, center, n, r200):
        in_halo = nbe.all_inside(self.pos3d, center, n*r200)
        self.pos3d = self.pos3d[in_halo] - center
        self.mass = self.mass[in_halo]
        self.pot = self.pot[in_halo]
        self.hsml = self.hsml[in_halo]
        self.vel3d = self.vel3d[in_halo]
        self.id = self.id[in_halo]
        self.R = np.sqrt((self.pos3d[:,0]**2)+(self.pos3d[:,1]**2))
        self.r = np.sqrt((self.pos3d[:,0]**2)+(self.pos3d[:,1]**2)+(self.pos3d[:,2]**2))
        ### velocities ###
        self.v = np.sqrt((self.vel3d[:,0]**2) + (self.vel3d[:,1]**2) + (self.vel3d[:,2]**2))
        self.vR = (self.vel3d[:,0]*self.pos3d[:,0] + self.vel3d[:,1]*self.pos3d[:,1])/ self.R
        self.vr = (self.vel3d[:,0]*self.pos3d[:,0] + self.vel3d[:,1]*self.pos3d[:,1] + self.vel3d[:,2]*self.pos3d[:,2])/ self.r
        self.vphi = (-self.vel3d[:,0]*self.pos3d[:,1] + self.vel3d[:,1]*self.pos3d[:,0] )/ self.R
        #### other params ###
        self.total_m =  np.sum(self.mass)

    def rotate(self,T):
        pos = self.pos3d
        self.pos3d = nbe.matrix_vs_vector(T,pos)
        self.vel3d = nbe.matrix_vs_vector(T,self.vel3d)


    def shift(self,center):
        self.pos3d = self.pos3d - center
        self._center_history = np.vstack((self._center_history,center))

