import sys
import math
import glob
import cmath
#import cfalcon
import subprocess
import numpy as np
from unsio import *
import scipy.special as sp
from numpy import exp, sqrt
import nbody_essentials as nbe
import scipy.integrate as integrate
from sklearn.neighbors import KDTree
#CF =cfalcon.CFalcon()

class Clumps: 
    def __init__(self, file_path,p, **kwargs):
        self.p = p 
        comov= kwargs.get('comov',False)
        self._center_history = np.array([[0,0,0]])
        self.data = nbe._get_center(file_path, clumps=True)
        if (comov):
            self.pos3d = self.data[:,4:7] * p.simutokpc / self.p.aexp
        else:
            self.pos3d = self.data[:,4:7] * p.simutokpc

        self.mass = self.data[:,10] * p.simutoMsun
        self.relevance = self.data[:,11] 
        self.n_cell = self.data[:,3] 
        self.level = self.data[:,1] 


    def halo_Only(self, center, n, r200):
        self.r = np.sqrt((self.pos3d[:,0]**2)+(self.pos3d[:,1]**2)+(self.pos3d[:,2]**2))
        in_halo = np.where(self.r <= n*r200)
        self.pos3d = self.pos3d[in_halo]
        self.mass = self.mass[in_halo]
        self.n_cell = self.n_cell[in_halo]
        self.relevance = self.relevance[in_halo]
        self.level = self.level[in_halo]
        self.data = self.data[in_halo]
        self.r = self.r[in_halo]
    
    def rotate(self,T):
        pos = self.pos3d
        self.pos3d = nbe.matrix_vs_vector(T,pos) 
    
    def shift(self, center):
        self.pos3d = self.pos3d - center
        self._center_history = np.vstack((self._center_history,center))
