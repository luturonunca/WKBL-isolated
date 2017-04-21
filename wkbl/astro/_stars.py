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
         
        if self.uns.isValid()!=True:
            sys.exit("\n\n\n\n\n\n\nERROR:"+file_path+" is not a valid file !!!!!\n\n\n\n\n")
        ok = self.uns.nextFrame("")
        ok, pos = self.uns.getArrayF("all","pos")
        ok, self.age = self.uns.getArrayF("stars","age")
        ok, vel = self.uns.getArrayF("all","vel")
        ok, mass = self.uns.getArrayF("all","mass")
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
        print len(self.pos3d)
        self.vel3d = vel.reshape(len(vel)/3,3)
        self.mass = mass * self._p.simutoMsun
        ## centering to the stars center ##
        try:
            if len(self.pos3d):
                tree = KDTree(np.squeeze(self.pos3d))
                in_halo = tree.query_radius(center,r_search)[0]
                pos_around_phew = self.pos3d[in_halo]
                m_around_phew = self.mass[in_halo]
                self.center_com = nbe.real_center(pos_around_phew,m_around_phew)
            else:
                self.center_com = np.zeros(3)
        except:
            sys.exit("ERROR in tree")
 

    def halo_Only(self, center, n, r200):
        in_halo = nbe.all_inside(self.pos3d, center,n*r200)
        self.pos3d = self.pos3d[in_halo] - center
        self.mass = self.mass[in_halo]
        self.age = self.age[in_halo]
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
        pos = self.pos3d
        self.pos3d = nbe.matrix_vs_vector(T,pos)

    def shift(self, center):
        self.pos3d = self.pos3d - center

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

       
