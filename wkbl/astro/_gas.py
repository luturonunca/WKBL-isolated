from . import component as comp
import numpy as np

class _gas(comp.Component):
    def __init__(self, file_path,p,dens=True, **kwargs):
        self.get_sigma = kwargs.get('virial',p.nmlexist)
        super().__init__(file_path,"gas",p)
        ok, temp = self.uns.getData("gas","temp")
        ok, rho = self.uns.getData("gas","rho")
        ok, self.rho1 = self.uns.getData("gas","rho")
        ok, self.pot = self.uns.getData("gas","pot")
        try:
            ok, self.B_left_x  = self.uns.getData("hydro","4")
            ok, self.B_left_y  = self.uns.getData("hydro","5")
            ok, self.B_left_z  = self.uns.getData("hydro","6")
            ok, self.B_right_x = self.uns.getData("hydro","7")
            ok, self.B_right_y = self.uns.getData("hydro","8")
            ok, self.B_right_z = self.uns.getData("hydro","9")
            self.bx = 0.5*(self.B_left_x+self.B_right_x)
            self.by = 0.5*(self.B_left_y+self.B_right_y)
            self.bz = 0.5*(self.B_left_z+self.B_right_z)
            self.bnorm = np.sqrt(self.bx**2 + self.by**2 + self.bz**2)
            self.va  = self.bnorm/rho*self._p.simutokms
            ok, self.non_th_pres = self.uns.getData("hydro","10")#*self._p.simutoErgscm3
            ok, self.pres = self.uns.getData("hydro","11")#*self._p.simutoErgscm3
            temp2 = self.pres/rho
        except:
            flag=1
        self.rho =  rho * self._p.unitd*(self._p.kpctocm**3)/1e3/self._p.msuntokg
        self.non_th_pres *= self._p.simutoErgscm3
        self.pres *= self._p.simutoErgscm3
        ok, hsml = self.uns.getData("gas","hsml")
        self.hsml = hsml * self._p.simutokpc
        shift1 = (np.random.rand(len(hsml))-0.5)*self.hsml
        shift2 = (np.random.rand(len(hsml))-0.5)*self.hsml
        shift3 = (np.random.rand(len(hsml))-0.5)*self.hsml
        self.pos3d[:,0]+=shift1
        self.pos3d[:,1]+=shift2
        self.pos3d[:,2]+=shift3
        self.tokelvin = self._p.mH / (1.3806200e-16) * (self._p.unitl / self._p.unitt)**2
        self.temp = temp * self.tokelvin
        if (self.get_sigma):
            ok, sigma = self.uns.getData("hydro",str(self._p.nener))
            self.sigma2 = sigma*(self._p.simutokms**2)
            #self.cs2 = (1.6667-1.) * self.pres * (self._p.simutokms**2)
            g_star = 1.6
            
    def halo_Only(self, center,n , r200,simple=False):
        super().halo_Only(center,n , r200,simple=simple)
        in_halo = np.where(self.r <= n*r200)
        self.pot = self.pot[in_halo]
        self.temp = self.temp[in_halo]
        self.pres = self.pres[in_halo]
        self.non_th_pres = self.non_th_pres[in_halo]
        self.hsml = self.hsml[in_halo]
        self.rho = self.rho[in_halo]
        self.B_left_x = self.B_left_x[in_halo]
        self.B_left_y = self.B_left_y[in_halo]
        self.B_left_z = self.B_left_z[in_halo]
        self.B_right_x = self.B_right_x[in_halo]
        self.B_right_y = self.B_right_y[in_halo]
        self.B_right_z = self.B_right_z[in_halo]
        self.bx = self.bx[in_halo]
        self.by = self.by[in_halo]
        self.bz = self.bz[in_halo]
        self.bR = (self.bx*self.pos3d[:,0] + self.by*self.pos3d[:,1])/ self.R
        self.bphi = (-self.bx*self.pos3d[:,1] + self.by*self.pos3d[:,0] )/ self.R
        if (self.get_sigma):
            self.met = self.met[in_halo]
            self.sigma2 = self.sigma2[in_halo]
            self.cs2 = self.cs2[in_halo]
        self.r = self.r[in_halo]
        
