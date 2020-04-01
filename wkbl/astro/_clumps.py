import numpy as np
from . import nbody_essentials as nbe
import glob

class Clumps: 
    def __init__(self, file_path,p, **kwargs):
        self.p = p 
        comov= kwargs.get('comov',False)
        rs= kwargs.get('rs',"")
        if len(rs)>1:self.ROCKSTAR=True
        self._center_history = np.array([[0,0,0]])
        self.data = nbe._read_extra(file_path,path=rs, clumps=True,rockstar=self.ROCKSTAR)
        if (self.ROCKSTAR):
            h = p.h 
            files = glob.glob(file_path+"/halos*.ascii")
            x = y = z = vx = vy = vz = rvir = m = np.array([])
            self.pos3d = np.zeros((len(self.data[:,0]),3))
            self.vel3d = np.zeros((len(self.data[:,0]),3))
            self.pos3d[:,0] = self.data[:,8]*1e3/h 
            self.pos3d[:,1] = self.data[:,9]*1e3/h
            self.pos3d[:,2] = self.data[:,10]*1e3/h
            self.vel3d[:,0] = self.data[:,11]/h
            self.vel3d[:,1] = self.data[:,12]/h
            self.vel3d[:,2] = self.data[:,13]/h
            self.mass = self.data[:,2]/h**2
            self.rvir = self.data[:,4]/h
            if not (comov):
                self.pos3d = self.pos3d * p.aexp
        else:
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
        if (self.ROCKSTAR):
            self.vel3d = self.vel3d[in_halo]
            self.rvir = self.rvir[in_halo]
        else:
            self.n_cell = self.n_cell[in_halo]
            self.relevance = self.relevance[in_halo]
            self.level = self.level[in_halo]
            self.data = self.data[in_halo]
        self.r = self.r[in_halo]
    
    def rotate(self,T):
        pos = self.pos3d
        self.pos3d = nbe.matrix_vs_vector(T,pos) 
        if (self.ROCKSTAR):
            vel = self.vel3d
            self.vel3d = nbe.matrix_vs_vector(T,vel) 

    
    def shift(self, center):
        self.pos3d = self.pos3d - center
        self._center_history = np.vstack((self._center_history,center))
