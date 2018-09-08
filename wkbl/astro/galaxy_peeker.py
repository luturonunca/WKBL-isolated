import sys
import math
import glob
import cmath
import subprocess
import numpy as np
from py_unsio import *
import scipy.special as sp
from numpy import exp, sqrt
import scipy.integrate as integrate
from sklearn.neighbors import KDTree
import nbody_essentials as nbe
from _dark_matter import _dark_matter
from _stars import _stars
from _gas import _gas
import datetime
import warnings
warnings.filterwarnings('ignore')

############### Unified  ################# 


class Galaxy_Hound:
    def __init__(self, file_path,getcen=False,**kwargs):
        # get them vars !!!!ONLy RAMSES FOR NOW
        self.file = file_path
        self.dmo = False
        halonu = len(glob.glob(file_path+"/halo*"))
        if halonu==0:
            newage = False
        else:
            newage = True
        gas = kwargs.get('gas',True)
        self.quiet = kwargs.get('quiet',False)
        self.p = nbe.Info_sniffer(file_path, newage=newage)
        hsml = kwargs.get('hsml',False)
        self.flush = kwargs.get('flush',False)
        dens = kwargs.get('dens',False)
        comov = kwargs.get('comov',True)
        self._dms, self._sts, self._gss  = False, False, False
        halo_vel = kwargs.get('halo_vel',[0.,0.,0.])    ##########
        #loadDatat
        self.n_tot, self.n_dm, self.n_st = nbe.check_particles(file_path)
        if self.n_dm > 0:
            if not self.quiet: print "loading Dark matter.."
            self.dm = _dark_matter(file_path, self.p,comov=comov)
            self._dms = True
        if self.n_st > 0:
            if not self.quiet: print "loading Stars.."
            self.st = _stars(file_path, self.p, comov=comov)
            self._sts = True
            if  gas==True:    
                if not self.quiet: print "loading Gas.."
                self.gs = _gas(file_path, self.p,comov=comov)
                self._gss = True
        else:
            self.dmo = True

        if (getcen):
            po, ma = np.copy(self.st.pos3d), np.copy(self.st.mass)
            centro_com_st = nbe.real_center(po,ma,n=7000)
            self.center_shift(centro_com_st)

    def r_virial(self, r_max,r_min=0,rotate=True,n=2.5,bins=512,quiet=False):
        positions = np.array([], dtype=np.int64).reshape(0,3)
        masses = np.array([], dtype=np.int64)
        if (self._dms):
            positions = np.vstack([positions,self.dm.pos3d])
            masses = np.append(masses,self.dm.mass)
        if (self._sts):
            positions = np.vstack([positions,self.st.pos3d])
            masses = np.append(masses,self.st.mass)
        if (self._gss):
            positions = np.vstack([positions,self.gs.pos3d])
            masses = np.append(masses,self.gs.mass)
        
        r = np.sqrt((positions[:,0])**2 +(positions[:,1])**2 +(positions[:,2])**2 )
        mhist, rhist = np.histogram(r,range=(0.0,r_max),bins=bins, weights=masses )
        vol_bin = (4./3.)*np.pi*(rhist[:-1]**3)
        r_bin = rhist[:-1]+ 0.5*(rhist[2]-rhist[1])
        rho_s = np.cumsum(mhist) / vol_bin
        self.r200 = r_bin[np.argmin(np.abs(rho_s - (200 * self.p.rho_crit)))]
        self.r97 = r_bin[np.argmin(np.abs(rho_s - (97 * self.p.rho_crit)))]
        rnot = False

        if (rotate)and(self._sts):
            if (self.flush):self.redefine(n,simple=True)
            if self.p.Z>2:
                self.rotate_galaxy(rmin=0.5)
            else:
                self.rotate_galaxy()
            self.redefine(n)
            D = np.dot(self.matrix_T,np.dot(self.matrix_P,np.transpose(self.matrix_T)))
            if not self.quiet:
                print '| r_200 = {0:.2f}'.format(self.r200)
                print '| Diagonal matrix computed '
                print '|    | {0}, {1}, {2}|'.format(int(D[0,0]),int(D[0,1]),int(D[0,2]))
                print '| D =| {0}, {1}, {2}|'.format(int(D[1,0]),int(D[1,1]),int(D[1,2]))
                print '|    | {0},  {1}, {2}|'.format(int(D[2,0]),int(D[2,1]),int(D[2,2]))
        elif (self.dmo): 
            self.redefine(n)
    
    def center_shift(self,nucenter):
        self.center = np.zeros(3)
        if (self._dms):
            self.dm.shift(nucenter)
        if (self._sts):
            self.st.shift(nucenter)
        if (self._gss):
            self.gs.shift(nucenter)
      
    def redefine(self,n,simple=False):
        if (self._dms):
            self.dm.halo_Only(self.center, n, self.r200,simple=simple)
        if (self._sts):
            self.st.halo_Only(self.center, n, self.r200, self.r97,simple=simple)
        if (self._gss):
            self.gs.halo_Only(self.center, n, self.r200,simple=simple)

    def rotate_galaxy(self,rmin=3,rmax=10):
        r2 = (self.st.pos3d[:,0])**2 +(self.st.pos3d[:,1])**2 +(self.st.pos3d[:,2])**2
        pos_ring = self.st.pos3d[(r2<rmax**2)&(r2>rmin**2)]
        P = np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                first = np.average(pos_ring[:,i]*pos_ring[:,j])
                second =(np.average(pos_ring[:,i])*np.average(pos_ring[:,j]))
                P[i][j] = first - second
        eigen_values,evecs = np.linalg.eig(P)
        order = np.argsort(abs(eigen_values))
        T = np.zeros((3,3))
        T[0],T[1],T[2] = evecs[:,order[2]],evecs[:,order[1]],evecs[:,order[0]]
        self.matrix_T = T
        self.matrix_P = P
        if (self._dms):
            self.dm.rotate(T)
        if (self._sts):        
            self.st.rotate(T)        
        if (self._gss):
            self.gs.rotate(T)
                
   
    def save_galaxy(self, name, fltype, component):
        unsout=CunsOut(name,fltype)
        ages = False
        if component == "stars":
            length = len(self.st.pos3d)
            pos_out = self.st.pos3d.reshape(length*3).astype(np.float32, copy=False)
            mass_out = self.st.mass.astype(np.float32, copy=False)
            age_out = self.st.age.astype(np.float32, copy=False)
            ages = True
        if component == "halo":
            length = len(self.dm.pos3d)
            pos_out = self.dm.pos3d.reshape(length*3).astype(np.float32, copy=False)
            mass_out = self.dm.mass.astype(np.float32, copy=False)
        if component == "gas":
            length = len(self.gs.pos3d)
            pos_out = self.gs.pos3d.reshape(length*3).astype(np.float32, copy=False)
            mass_out = self.gs.mass.astype(np.float32, copy=False)
        ok=unsout.setArrayF(component,"pos",pos_out)
        ok=unsout.setArrayF(component,"mass",mass_out)
        if (ages):
            ok=unsout.setArrayF(component,"age",age_out)
        unsout.save()
   
    def get_shell(self,lim_min,lim_max):
        shell_cond = (self.r<=lim_max)&(self.r>=lim_min)
        self.v_shell = self.v[shell_cond]
        self.vr_shell = self.vr[shell_cond]
        self.vR_shell = self.vR[shell_cond]
        self.vphi_shell = self.vphi[shell_cond]
        self.vz_shell = self.vz[shell_cond]

    def get_ring_xy(self,lim_min,lim_max, thickness):
        ringXY_cond = (self.r<=lim_max)&(self.r>=lim_min)&(abs(self.z)<=thickness)
        self.v_ringXY = self.v[ringXY_cond]
        self.vr_ringXY = self.vr[ringXY_cond]
        self.vR_ringXY = self.vR[ringXY_cond]
        self.vphi_ringXY = self.vphi[ringXY_cond]
        self.vz_ringXY = self.vz[ringXY_cond]
        self.x_ringXY = self.x[ringXY_cond]
        self.y_ringXY = self.y[ringXY_cond]
        self.z_ringXY = self.z[ringXY_cond]
        self.vx_ringXY = self.vx[ringXY_cond]
        self.vy_ringXY = self.vy[ringXY_cond]

    def get_ring_yz(self,lim_min,lim_max, thickness):
        ringYZ_cond = (self.r<=lim_max)&(self.r>=lim_min)&(abs(self.x)<=thickness)
        self.v_ringYZ = self.v[ringYZ_cond]
        self.vr_ringYZ = self.vr[ringYZ_cond]
        self.vR_ringYZ = self.vR[ringYZ_cond]
        self.vphi_ringYZ = self.vphi[ringYZ_cond]
        self.vz_ringYZ = self.vz[ringYZ_cond]

    def get_ring_zx(self,lim_min,lim_max, thickness):
        ringZX_cond = (self.r<=lim_max)&(self.r>=lim_min)&(abs(self.y)<=thickness)
        self.v_ringZX = self.v[ringZX_cond]
        self.vr_ringZX = self.vr[ringZX_cond]
        self.vR_ringZX = self.vR[ringZX_cond]
        self.vphi_ringZX = self.vphi[ringZX_cond]
        self.vz_ringZX = self.vz[ringZX_cond]

    def get_dens(self,lim_min,lim_max, thickness,plane="x"):
        shell_cond = (self.r<=lim_max)&(self.r>=lim_min)
        if plane == "x":
            plane_cond = (abs(self.z)<=thickness)
        elif plane == "y":
            plane_cond = (abs(self.x)<=thickness)
        elif plane == "z":
            plane_cond = (abs(self.y)<=thickness)
        else:
            plane_cond = (True)
        condition = shell_cond & plane_cond
        self.rho_ring = self.rho[condition]
        self.rho_shell = self.rho[shell_cond]
