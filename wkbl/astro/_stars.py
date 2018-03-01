import sys
import math
import glob
import cmath
import subprocess
import numpy as np
from py_unsio import *
import nbody_essentials as nbe
import scipy.special as sp
from numpy import exp, sqrt
import scipy.integrate as integrate
from sklearn.neighbors import KDTree
from _sf_info import SF_info

class _stars:
    def __init__(self, file_path,p, **kwargs):
        self._p = p
        center = kwargs.get('center',[0.,0.,0.])    ##########
        self.uns = CunsIn(file_path,"stars","all",False)
        hsml = kwargs.get('hsml',False)
        dens = kwargs.get('dens',False)
        comov = kwargs.get('comov',False)
        r_search = kwargs.get('r_search',200.)
        self.halo_vel = kwargs.get('halo_vel',[0.,0.,0.])    ##########
        try:
            self.sf_info = SF_info(file_path,p,comov=comov)
            self.gotsfInfo = True
        except:
            self.gotsfInfo = False
 
        if self.uns.isValid()!=True:
            sys.exit("\n\n\n\n\n\n\nERROR:"+file_path+" is not a valid file !!!!!\n\n\n\n\n")
        ok = self.uns.nextFrame("")
        ok, pos = self.uns.getArrayF("all","pos")
        ok, age = self.uns.getArrayF("stars","age")
        ok, vel = self.uns.getArrayF("all","vel")
        ok, mass = self.uns.getArrayF("all","mass")
        ok, self.metal = self.uns.getArrayF("all","metal")
        ok, self.id = self.uns.getArrayI("all","id")
        if dens ==True:
            ok, rho = self.uns.getArrayF("all","rho")
            self.rho =  rho * simutoGeVcm3
        if hsml ==True:
            ok, self.hsml = self.uns.getArrayF("all","hsml")
        ### coordinates ###
        pos = pos * self._p.simutokpc
        vel = vel * self._p.simutokms
        if (comov):
            self.pos3d = pos.reshape(len(pos)/3,3) / self._p.aexp
        else:
            self.pos3d = pos.reshape(len(pos)/3,3)
        self.vel3d = vel.reshape(len(vel)/3,3)
        self.mass = mass * self._p.simutoMsun
        self.age = age *self._p.unitt / (3600.*24.*365*1e9) # stars age to Gyrs
 

    def halo_Only(self, center, n, r200,r97):
        #### sf history ###
        if (self.gotsfInfo):
            self.sf_info.halo_Only(center, n, r200)
        in_halo = nbe.all_inside(self.pos3d, center,n*r200)
        self.pos3d = self.pos3d[in_halo] - center
        self.mass = self.mass[in_halo]
	self.age = self.age[in_halo]
        self.metal = self.metal[in_halo]
        self.vel3d = self.vel3d[in_halo]
        self.id = self.id[in_halo]
        self.R = np.sqrt((self.pos3d[:,0]**2)+(self.pos3d[:,1]**2))
        self.r = np.sqrt((self.pos3d[:,0]**2)+(self.pos3d[:,1]**2)+(self.pos3d[:,2]**2))
        ### unitary directions###
        po, R = self.pos3d, np.column_stack([self.R,self.R,self.R])
        ur = po / np.column_stack([self.r,self.r,self.r])
        uR = po / R
        uphi = (po * [-1,1,1]) / np.column_stack([self.r,self.r,self.r]) 
        thetax = (uphi[:,1]*ur[:,2])-(uphi[:,2]*ur[:,1])
        thetay = (uphi[:,2]*ur[:,0])-(uphi[:,0]*ur[:,2])
        thetaz = (uphi[:,0]*ur[:,1])-(uphi[:,1]*ur[:,0])
        utheta = np.column_stack([thetax,thetay,thetaz])
        ### velocities ###
        average_v = np.array([np.mean(self.vel3d[:,0]),np.mean(self.vel3d[:,1]),np.mean(self.vel3d[:,2])])
        self.vel3d = self.vel3d - average_v
        vx,vy,vz = self.vel3d[:,0],self.vel3d[:,1],self.vel3d[:,2]
        self.v = np.sqrt((vx**2) + (vy**2) + (vz**2))

        self.vR = self.vel3d * uR#(vx*self.pos3d[:,0] + vy*self.pos3d[:,1])/ self.R
        self.vr = self.vel3d * ur#(vx*self.pos3d[:,0] + vy*self.pos3d[:,1] + vz*self.pos3d[:,2])/ self.r
        self.vphi = np.sum(self.vel3d * uphi, axis=1)#(-vx*self.pos3d[:,1] + vy*self.pos3d[:,0] )/ self.R
        self.vtheta = np.sum(self.vel3d * utheta, axis=1)#(self.vR*self.pos3d[:,2] - vz*self.R) / self.r
        #### other params ###
        self.total_m =  np.sum(self.mass)
        self.fire_m , self.fire_r= nbe.FIRE_st_mass(self.mass,self.r,r97)

    def Age_cut(self, age_cut):
        in_gal = (self.age < age_cut)
        self.pos3d = self.pos3d[in_gal] 
        self.mass = self.mass[in_gal]
        self.vel3d = self.vel3d[in_gal]
        self.id = self.id[in_gal]
        self.R = self.R[in_gal]
        self.r = self.r[in_gal]
        ### velocities ###
        self.v = self.v[in_gal]
        self.vR = self.vR[in_gal]
        self.vr = self.vr[in_gal]
        self.vphi = self.vphi[in_gal]
        #### other params ###
        self.total_m =  np.sum(self.mass)
        self.age = self.age[in_gal]
        
    def rotate(self,T):
        if (self.gotsfInfo):
            self.sf_info.rotate(T)
        pos = self.pos3d
        self.pos3d = nbe.matrix_vs_vector(T,pos)

    def shift(self, center):
        self.pos3d = self.pos3d - center
        if (self.gotsfInfo):
            self.sf_info.shift(center)

    def get_M_virS(self,r200,r97):
	# total mass inside r200
        print " M_{r_200}.."
        self.M200 = np.sum(self.mass[(self.r<=r200)])
	# total mass inside r97
        print " M_{r_97}.."
        self.M97 = np.sum(self.mass[(self.r<=r97)])
        # setellar mass as defined in Fire 2 [arxiv:1702.06148v1] (p.16 footnote 9)
        print " M_{Fire2}.."
        r = self.r
        mass_1 = self.mass[(r < 0.15 * r97)]
        r_1 = self.r[(r < 0.15 * r97)]
	r_half_1 = nbe.half_mass(mass_1,r_1)
	mass_2 = self.mass[(r<3*r_half_1)]
	r_2 = self.r[(r<3*r_half_1)]
	self.r_F2 = nbe.half_mass(mass_2, r_2)
        self.M_F2 = np.sum(self.mass[(self.r<self.r_F2)]) 

       
