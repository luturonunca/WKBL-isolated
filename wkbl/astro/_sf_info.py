import sys,os
import math
import glob
import cmath
#import cfalcon
import subprocess
import numpy as np
from unsio import *
import scipy.special as sp
from numpy import exp, sqrt
from . import nbody_essentials as nbe
import scipy.integrate as integrate
from sklearn.neighbors import KDTree
#CF =cfalcon.CFalcon()

class SF_info: 
    def __init__(self, file_path,p, **kwargs):
        # load st_info where the current stars where formed
        star_files = glob.glob(file_path+"/star*")
        if len(star_files)==0:
            # if there is no star files in the output directory
            # proceed to use the previous output, sometimes the 
            # very last output does not present any star files because
            # is just a final closing output
            num = float(file_path[-5:])
            if num<100:
                # load de stars info from the previous output
                file_path = file_path[:-2]+str(int(file_path[-2:])-1)
            else:
                file_path = file_path[:-3]+str(int(file_path[-3:])-1)
        self.p = p 
        comov= kwargs.get('comov',False)
        self._center_history = np.array([[0,0,0]])
        self.data = nbe._read_extra(file_path, sf_hist=True)
        stars = np.where(self.data[:,0]==0)[0]
        if (comov):
            self.pos3d = self.data[stars,4:7] * p.simutokpc / p.aexp
        else:
            self.pos3d = self.data[stars,4:7] * p.simutokpc
        self.mass_st = self.data[stars,3] * p.simutoMsun
        self.id = self.data[stars,1]
        self.hsml = 25000./(2.**(self.data[stars,2]))
        rho = self.data[stars,10]
        self.rho = rho*self.p.simutoMsun / (self.p.simutokpc**3)
        self.tokelvin = 1.66e-27 / (1.3806200e-19) * (self.p.unitl / self.p.unitt)**2
        self.B_left_x = self.data[:,14]
        self.B_left_y = self.data[:,15]
        self.B_left_z = self.data[:,16]
        self.B_right_x = self.data[:,17]
        self.B_right_y = self.data[:,18]
        self.B_right_z = self.data[:,19]
        self.bx = 0.5*(self.B_left_x+self.B_right_x)
        self.by = 0.5*(self.B_left_y+self.B_right_y)
        self.bz = 0.5*(self.B_left_z+self.B_right_z)
        self.bnorm = np.sqrt(self.bx**2 + self.by**2 + self.bz**2)
        self.va  = self.bnorm/rho*self.p.simutokms
        self.non_th_pres = self.data[stars,20]*self.p.simutoErgscm3
        self.pres = self.data[stars,21]*self.p.simutoErgscm3
        self.temp = (self.data[stars,21]/rho)*self.tokelvin
        self.size = self.p.boxlen*1e5/2**self.data[stars,2] 
        self.level = self.data[stars,2] 
        self.mass = self.rho*(self.size**3)
         
    def halo_Only(self, center, n, r):
        self.r = np.sqrt((self.pos3d[:,0]**2)+(self.pos3d[:,1]**2)+(self.pos3d[:,2]**2))
        in_halo = np.where(self.r <= n*r)
        self.pos3d = self.pos3d[in_halo] - center
        self.id = self.id[in_halo]
        self.rho = self.rho[in_halo]
        self.size = self.size[in_halo]
        self.level = self.level[in_halo]
        self.temp = self.temp[in_halo]
        self.data = self.data[in_halo]
        self.r = self.r[in_halo]

    
    def rotate(self,T):
        pos = self.pos3d
        self.pos3d = nbe.matrix_vs_vector(T,pos) 
   
    def id_sort(self):
        sorted_by_id = np.argsort(self.id)
        self.pos3d = self.pos3d[sorted_by_id]
        self.id = self.id[sorted_by_id]
        self.rho = self.rho[sorted_by_id]
        self.rho_crit = self.rho_crit[sorted_by_id]
        self.temp2 = self.temp2[sorted_by_id]
        self.met = self.met[sorted_by_id]
        self.alpha0= self.alpha0[sorted_by_id]
        self.sigma2 = self.sigma2[sorted_by_id]
        self.cs2 = self.cs2[sorted_by_id]
        self.M2 = self.M2[sorted_by_id]    
    
  
    def shift(self, center):
        self.pos3d = self.pos3d - center
        self._center_history = np.vstack((self._center_history,center))

