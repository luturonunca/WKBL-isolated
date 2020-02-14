from . import Component
import numpy as np

class _stars(Component):
    def __init__(self, file_path,p, **kwargs):
        dens = kwargs.get('dens',False)
        comov = kwargs.get('comov',False)
        r_search = kwargs.get('r_search',200.)
        try:
            self.sf_info = SF_info(file_path,p,comov=comov)
            self.gotsfInfo = True
        except:
            self.gotsfInfo = False
        super().__init__(file_path,"stars",p)
        ok, age = self.uns.getArrayF("stars","age")
        ok, self.metal = self.uns.getArrayF("stars","metal")
        ok, self.id = self.uns.getArrayI("all","id")
        self.age = age *self._p.unitt / (3600.*24.*365*1e9) / self._p.aexp**2 # stars age to Gyrs
        
    def halo_Only(self, center, n, r200, simple=False):
        #### sf history ####
        if (self.gotsfInfo):
            self.sf_info.halo_Only(center, n, r200)
        super().halo_Only(center,n , r200, simple=simple)
        in_halo = np.where(self.r <= n*r200)
        self.age = self.age[in_halo]
        self.metal = self.metal[in_halo]
        self.r = self.r[in_halo]
 
    def shift(self,center):
        if (self.gotsfInfo):self.sf_info.shift(center)
        super().shift(center)
        
