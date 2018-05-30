import sys
import math
import glob
import cmath
#import cfalcon
import subprocess
import numpy as np
from py_unsio import *
import scipy.special as sp
from numpy import exp, sqrt
import nbody_essentials as nbe
import scipy.integrate as integrate
#CF =cfalcon.CFalcon()
#from _clumps import Clumps


class _dark_matter:
    def __init__(self, file_path,p, **kwargs):
        self._p = p
        self.file = file_path
        self._dens = False
        self._center_history = np.array([0.,0.,0.])    ##########
        self.uns = CunsIn(file_path,"halo","all",False)
        hsml = kwargs.get('hsml',False)
        dens = kwargs.get('dens',False)
        comov = kwargs.get('comov',False)
        #self.Clumps = Clumps(file_path,p,comov=comov) 
        self.halo_vel = kwargs.get('halo_vel',[0.,0.,0.])    ##########

        if self.uns.isValid()!=True:
            sys.exit("\n\n\n\n\n\n\nERROR:"+file_path+" is not a valid file !!!!!\n\n\n\n\n")
        ok = self.uns.nextFrame("")
        ok, pos = self.uns.getArrayF("all","pos")
        ok, vel = self.uns.getArrayF("all","vel")
        ok, mass = self.uns.getArrayF("all","mass")
        ok, self.id = self.uns.getArrayI("all","id")
        if dens ==True:
            ok, rho = self.uns.getArrayF("all","rho")
            self.rho =  rho * simutoGeVcm3
        ### coordinates ###
        if (comov):
            pos = pos * self._p.simutokpc / self._p.aexp
        else:
            pos = pos * self._p.simutokpc
        vel = vel * self._p.simutokms
        self.pos3d = pos.reshape(len(pos)/3,3)
        self.vel3d = vel.reshape(len(vel)/3,3)
        self.mass = mass * self._p.simutoMsun
         
    def halo_Only(self, center,n , r200):
        #### clumps ###
        #self.Clumps.halo_Only(center, n, r200)
        ### dm particles ##
        self.r = np.sqrt((self.pos3d[:,0]**2)+(self.pos3d[:,1]**2)+(self.pos3d[:,2]**2))
        in_halo = np.where(self.r <= n*r200)
        self.pos3d = self.pos3d[in_halo]
        self.mass = self.mass[in_halo]
        # substract average velocities
        self.vel3d = self.vel3d[in_halo]
        average_v = np.array([np.mean(self.vel3d[:,0]),np.mean(self.vel3d[:,1]),np.mean(self.vel3d[:,2])])
        self.vel3d = self.vel3d - average_v
        self.id = self.id[in_halo]
        # spherical/cylindrical coordinates
        self.R = np.sqrt((self.pos3d[:,0]**2)+(self.pos3d[:,1]**2))
        self.r = self.r[in_halo]
        self.phi = np.arctan2(np.copy(self.pos3d[:,1]),np.copy(self.pos3d[:,0]))
        self.theta = np.arccos(np.copy(self.pos3d[:,0]),np.copy(self.r))
        ### velocities ###
        vx,vy,vz = self.vel3d[:,0],self.vel3d[:,1],self.vel3d[:,2]
        self.v = np.sqrt((vx**2) + (vy**2) + (vz**2))
        self.vR = (vx*self.pos3d[:,0] + vy*self.pos3d[:,1])/ self.R
        self.vr = (vx*self.pos3d[:,0] + vy*self.pos3d[:,1] + vz*self.pos3d[:,2])/ self.r
        self.vphi = (-vx*self.pos3d[:,1] + vy*self.pos3d[:,0] )/ self.R
        self.vtheta = (self.vR*self.pos3d[:,2] - vz*self.R) / self.r
        #### other params ###
        self.total_m =  np.sum(self.mass)
        self.center_com = nbe.real_center(self.pos3d,self.mass)

    def rotate(self,T):
        #self.Clumps.rotate(T)
	pos = self.pos3d
        self.pos3d = nbe.matrix_vs_vector(T,pos)
        self.vel3d = nbe.matrix_vs_vector(T,self.vel3d)
        
    def get_shperical_coords(self):
        """
        gets theta and phi coordinates
        can only be run after the halo is defined
        i.e after halo_Only()
        """
        self.phi = np.arctan2(np.copy(self.pos3d[:,1]),np.copy(self.pos3d[:,0]))
        self.theta = np.arccos(np.copy(self.pos3d[:,0]),np.copy(self.r))



 
    def shift(self,center):
        self.pos3d = self.pos3d - center
        self._center_history = np.vstack((self._center_history,center))
        #self.Clumps.shift(center)    

    def density_profile(self, bins, limit):
        r_p = np.logspace(-0.5, np.log10(limit),bins)
        def sph_dens(r):
            """
            spherical density profile
            """
            total_mass = np.sum(self.mass[(self.r < r)])
            dens = 4 * total_mass / 3. / np.pi / r**3
            return dens

        get_shp_dens = np.vectorize(sph_dens)
        return r_p , get_shp_dens(r_p)

    def den(self, bin_num, r200,  quiet=False):
        """
        Density profile of DM halo
        returns two arrays:
        1- radius normalize to r200, bin_num long
        2- density of the shell that correspond to each radius bin, bin_num-1 long
        """
        if not quiet:
            print "| Beware of the center you are considering for the Dark Matter halo here,"
            print "| depending on the center Cusp could turn into Core or viceversa "
            print "| check nbody_essentials.real_center"
        if not self._dens:
            try:
                ok,self.rho,_= CF.getDensity(np.array(self.pos3d.reshape(len(self.pos3d)*3),dtype=np.float32), self.mass)
                self._dens = True
            except: 
                sys.exit("--density could not be computed--")

        bins = np.logspace(-0.6, np.log10(r200), bin_num)
        array = []
        for i in range(0,len(bins)-1):
            array.append(np.average(self.rho[(self.r>bins[i])&(self.r<bins[i+1])]))
        return bins/r200, array

    def get_M_virS(self,r200,r97):
        # total mass inside r200
        print " M_{r_200}.."
        self.M200 = np.sum(self.mass[(self.r<=r200)])
        # total mass inside r97
        print " M_{r_97}.."
        self.M97 = np.sum(self.mass[(self.r<=r97)])
