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
        self.cell = self.data[:,3]

    def halo_Only(self, center, n, r):
        in_halo = nbe.all_inside(self.pos3d, center, r)
        self.pos3d = self.pos3d[in_halo] - center
        self.mass = self.mass[in_halo]
    
    def rotate(self,T):
        pos = self.pos3d
        self.pos3d = nbe.matrix_vs_vector(T,pos) 
    
    def shift(self, center):
        self.pos3d = self.pos3d - center
        self._center_history = np.vstack((self._center_history,center))
