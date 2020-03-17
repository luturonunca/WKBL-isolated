<<<<<<< HEAD
"""
gas  and AMR routine
"""
import sys
import math
import glob
import cmath
import subprocess
=======
from . import component as comp
>>>>>>> python3Ver
import numpy as np

<<<<<<< HEAD
class _gas:
    def __init__(self, file_path,p, **kwargs):
        self._p = p
        self._center_history = np.array([0.,0.,0.]) 
        self.uns = CunsIn(file_path,"gas","all",False)
        dens = kwargs.get('dens',False)
        self.get_sigma = kwargs.get('virial',True)
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
            self.sigma2 = sigma*(self._p.simutokms**2) 
        
        ### coordinates ###
        vel = vel * self._p.simutokms
        if (comov):
            pos = pos * self._p.simutokpc / self._p.aexp
        else:
            pos = pos * self._p.simutokpc
        gamma = 1.666
        self.temp = (gamma-1.0) * temp * self._p.scale_T2 #Kelvin
        self.temp2 = temp2 * self._p.scale_T2# Kelvin/mu
        self.rho =  rho * self._p.scale_d / self._p.scale_d_gas
        self.pot = pot #* self._p.simutokms**2 
        self.pos3d = pos.reshape(len(pos)/3,3)
        self.center_rho_max = self.pos3d[np.where(self.rho==self.rho.max())][0]
        self.vel3d = vel.reshape(len(vel)/3,3)
        self.mass = mass * self._p.simutoMsun
        self.hsml = hsml * self._p.simutokpc

    def halo_Only(self, center, n, r200):
        self.r = np.sqrt((self.pos3d[:,0]**2)+(self.pos3d[:,1]**2)+(self.pos3d[:,2]**2))
        in_halo = np.where(self.r <= n*r200)
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
        self.r = self.r[in_halo] 
        self.phi = np.arctan2(np.copy(self.pos3d[:,1]),np.copy(self.pos3d[:,0]))
        self.theta = np.arccos(np.copy(self.pos3d[:,0]),np.copy(self.r))
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
        vel = self.vel3d
        self.pos3d = nbe.matrix_vs_vector(T,pos)
        self.vel3d = nbe.matrix_vs_vector(T,vel)

    def shift(self,center):
        self.pos3d = self.pos3d - center
        self._center_history = np.vstack((self._center_history,center))

=======
class _gas(comp.Component):
    def __init__(self, file_path,p,dens=True, **kwargs):
        self.get_sigma = kwargs.get('virial',p.nmlexist)
        super().__init__(file_path,"gas",p)
        ok, temp = self.uns.getData("gas","temp")
        ok, rho = self.uns.getData("gas","rho")
        ok, self.pot = self.uns.getData("gas","pot")
        ok, self.pres = self.uns.getData("hydro","4")
        ok, self.met = self.uns.getData("hydro","5")
        temp2 = self.pres/rho
        self.rho =  rho * self._p.simutoMsun / (self._p.simutokpc**3)

        ok, hsml = self.uns.getData("gas","hsml")
        self.hsml = hsml * self._p.simutokpc
        shift1 = (np.random.rand(len(hsml))-0.5)*self.hsml
        shift2 = (np.random.rand(len(hsml))-0.5)*self.hsml
        shift3 = (np.random.rand(len(hsml))-0.5)*self.hsml
        self.pos3d[:,0]+=shift1
        self.pos3d[:,1]+=shift2
        self.pos3d[:,2]+=shift3
        self.tokelvin = self._p.mH / (1.3806200e-16) * (self._p.unitl / self._p.unitt)**2
        self.temp = temp * self.tokelvin
        if (self.get_sigma):
            ok, sigma = self.uns.getData("hydro",str(self._p.nener))
            self.sigma2 = sigma*(self._p.simutokms**2)
            self.cs2 = (1.6667-1.) * self.pres * (self._p.simutokms**2)
            g_star = 1.6
            
    def halo_Only(self, center,n , r200,simple=False):
        super().halo_Only(center,n , r200,simple=simple)
        in_halo = np.where(self.r <= n*r200)
        self.met = self.met[in_halo]
        self.pot = self.pot[in_halo]
        self.temp = self.temp[in_halo]
        self.pres = self.pres[in_halo]
        self.hsml = self.hsml[in_halo]
        self.rho = self.rho[in_halo]
        if (self.get_sigma):
            self.sigma2 = self.sigma2[in_halo]
            self.cs2 = self.cs2[in_halo]
        self.r = self.r[in_halo]
        
>>>>>>> python3Ver
