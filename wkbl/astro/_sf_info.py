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
from sklearn.neighbors import KDTree
#CF =cfalcon.CFalcon()

class SF_info: 
    def __init__(self, file_path,p, **kwargs):
        # load st_info where the current stars where formed
        num = float(file_path[-5:])
        if num<100:
            file_path = file_path[:-2]+str(int(file_path[-2:])-1)
        else:
            file_path = file_path[:-3]+str(int(file_path[-3:])-1)
        self.p = p 
        comov= kwargs.get('comov',False)
        self._center_history = np.array([[0,0,0]])
        self.data = nbe._get_center(file_path, sf_hist=True)
        stars = np.where(self.data[:,0]==0)[0]
        if (comov):
            self.pos3d = self.data[stars,4:7] * p.simutokpc / p.aexp
        else:
            self.pos3d = self.data[stars,4:7] * p.simutokpc
        self.mass = self.data[stars,3] * p.simutoMsun
        self.id = self.data[stars,1]
        self.hsml = 25000./(2.**(self.data[stars,2]))
        self.rho = self.data[stars,10]*self.p.scale_d / self.p.scale_d_gas
        self.temp2 = self.data[stars,15]* self.p.scale_T2
        self.met = self.data[stars,16] 
        self.sigma2 = self.data[stars,17] * (self.p.simutokms**2)
        cs2 = 0.666*self.data[stars,15] * (self.p.simutokms**2)
        cs2[np.where(cs2<=0)] = 0.002
        self.cs2 = cs2
        self.M2 = self.sigma2 / self.cs2
        factG = 3. / 4. / 2. / np.pi * 0.3089 * self.p.aexp 
        self.alpha0 = 0.5 * self.sigma2 / np.pi / factG / self.rho / (self.hsml)**2
        phi_x = 0.19
        self.rho_crit = ((np.pi**2)/5.) * (phi_x**2) * self.alpha0 * self.M2
         
    def halo_Only(self, center, n, r):
        self.r = np.sqrt((self.pos3d[:,0]**2)+(self.pos3d[:,1]**2)+(self.pos3d[:,2]**2))
        in_halo = np.where(self.r <= n*r)
        self.pos3d = self.pos3d[in_halo] - center
        self.id = self.id[in_halo]
        self.rho = self.rho[in_halo]
        self.rho_crit = self.rho_crit[in_halo]
        self.temp2 = self.temp2[in_halo]
        self.met = self.met[in_halo]
        self.alpha0= self.alpha0[in_halo]
        self.sigma2 = self.sigma2[in_halo]
        self.cs2 = self.cs2[in_halo]
        self.M2 = self.M2[in_halo]
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

