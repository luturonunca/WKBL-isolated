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
        file_path = file_path[:-2]+str(int(file_path[-2:])-1)
        self.p = p 
        comov= kwargs.get('comov',False)
        self._center_history = np.array([[0,0,0]])
        self.data = nbe._get_center(file_path, sf_hist=True)
        if (comov):
            self.pos3d = self.data[:,4:7] * p.simutokpc / p.aexp
        else:
            self.pos3d = self.data[:,4:7] * p.simutokpc
        self.mass = self.data[:,3] * p.simutoMsun
        self.id = self.data[:,1]
        self.rho = self.data[:,10] * self.p.simutoMsun * (self.p.simutokpc)**-3
        self.temp = self.data[:,14]* p.scale_T2
        self.met = self.data[:,15] 
        self.sigmas = self.data[:,17] * p.simutokms 

    def halo_Only(self, center, n, r):
        in_halo = nbe.all_inside(self.pos3d, center, r)
        self.pos3d = self.pos3d[in_halo] - center
        self.id = self.id[in_halo]
        self.rho = self.rho[in_halo]
        self.temp = self.temp[in_halo]
        self.met = self.met[in_halo]
        self.sigmas = self.sigmas[in_halo]
    
    def rotate(self,T):
        pos = self.pos3d
        self.pos3d = nbe.matrix_vs_vector(T,pos) 
    
    def shift(self, center):
        self.pos3d = self.pos3d - center
        self._center_history = np.vstack((self._center_history,center))
