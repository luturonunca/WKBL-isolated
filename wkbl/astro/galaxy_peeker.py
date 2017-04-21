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

 
 


class Galaxy_Hound:
    def __init__(self, file_path, component,getcen=True,**kwargs):
        # get them vars !!!!ONLy RAMSES FOR NOW
        self.file = file_path
        self.p = nbe.Info_sniffer(file_path)
        hsml = kwargs.get('hsml',False)
        dens = kwargs.get('dens',False)
        comov = kwargs.get('comov',False)
        self._dms = False
        self._sts = False
        self._gss = False
        halo_vel = kwargs.get('halo_vel',[0.,0.,0.])    ##########
        #loadDatat
	if component=="all":
            sys.exit('we do not play with that card (component=all), be specific')
        self.comp = component.split(',')
        if "halo" in  component:
            self._dms = True
        if "gas" in  component:
            self._gss = True
        if "stars" in  component:
            self._sts = True

        if (self._dms):
            print "loading Dark matter.."
            self.dm = _dark_matter(file_path, self.p,comov=comov)
        if (self._gss):
            print "loading Gas.."
            self.gs = _gas(file_path, self.p,comov=comov)
        if (self._sts):
            print "loading Stars.."
            self.st = _stars(file_path, self.p, comov=comov)
            #self.st = _stars(file_path, self.p,center=self.center, comov=comov)
            #self.center_shift(self.st.center_com)
        #if (getcen):
        #    self.center = nbe._get_center(file_path) * self.p.simutokpc
        #else:
        #    self.center = np.zeros(3)
 
    def r_virial(self, r_max,r_min=0,rotate=True,n=2.5):
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
        try:
            tree = KDTree(np.squeeze(positions))
            in_halo = tree.query_radius(self.center,r_max)[0]
            pos_halo = positions[in_halo]
            m_in_halo = masses[in_halo]
        except:
            sys.exit("ERROR in tree")

        #r = np.sqrt((pos_halo[:,0]-self.center[0])**2 +(pos_halo[:,1]-self.center[1])**2 +(pos_halo[:,2]-self.center[2])**2 )
        r = np.sqrt((pos_halo[:,0])**2 +(pos_halo[:,1])**2 +(pos_halo[:,2])**2 )
        r = r[np.argsort(r)]
        mass_sorted = m_in_halo[np.argsort(r)]
        rho_local = 2*self.p.rho_crit * 97.
        i = np.where(r>r_min)[0][0]
        rnot = True #True until r_200 its founded
        if r_min==0:
            msu = 0
        else:
            msu = np.sum(mass_sorted[i-1])
        try:
            while rho_local >  self.p.rho_crit * 97.: 
                msu += mass_sorted[i]
                rho_local =  (3. /4. / np.pi) * msu / (r[i])**3
                i+=1
                if rho_local <= self.p.rho_crit * 200. and (rnot):
                    self.r200, rnot = r[i], False
 
        except:
            print "virial radius did not converged "
            sys.exit()
        self.r97 = r[i]
        if (rotate)and(self._sts):
            print '| r_200 = {0}'.format(self.r200)
            print '---- taking particles inside {0} * r200'.format(n)
            self.redefine(n)
            print '| number of praticles inside {0} * r200 '.format(n)
            if (self._dms):
                print '| dm mass       =  {0} M_sun'.format(self.dm.total_m)
                print '| p_dm_200      =  {0} particles'.format(len(self.dm.pos3d))
            if (self._sts):
                print '| stellar mass  =  {0} M_sun'.format(self.st.total_m)
                print '| p_st_200      =  {0} psrticles'.format(len(self.st.pos3d))
            if (self._gss):
                print '| gas mass      =  {0} M_sun'.format(self.gs.total_m)
                print '| p_gs_200      =  {0} particles'.format(len(self.gs.pos3d))
            print '---- rotating galaxy '
            self.rotate_galaxy()
            D = np.dot(self.matrix_T,np.dot(self.matrix_P,np.transpose(self.matrix_T)))
            print '| Diagonal matrix computed '
            print '|    |{0},{1},{2}|'.format(int(D[0,0]),int(D[0,1]),int(D[0,2]))
            print '| D =|{0},{1},{2}|'.format(int(D[1,0]),int(D[1,1]),int(D[1,2]))
            print '|    |{0},{1},{2}|'.format(int(D[2,0]),int(D[2,1]),int(D[2,2]))
        elif (rotate):
            self.redefine(n)
    
    def center_shift(self,nucenter):
        self.center = np.zeros(3)
        if (self._dms):
            self.dm.shift(nucenter)
        if (self._sts):
            self.st.shift(nucenter)
        if (self._gss):
            self.gs.shift(nucenter)
      
    def redefine(self,n):
        if (self._dms):
            print "%%flag dm"
            self.dm.halo_Only(self.center, n, self.r200)
        if (self._sts):
            print "%%flag st"
            self.st.halo_Only(self.center, n, self.r200)
        if (self._gss):
            print "%%flag gs"
            self.gs.halo_Only(self.center, n, self.r200)

    def rotate_galaxy(self,rmin=3,rmax=10):
        pos_ring = self.st.pos3d[(self.st.r<rmax)&(self.st.r>rmin)]
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
        #if((self._dms)and(self._sts)):
            #of_v = self.dm.center_com - self.st.center_com
            #self.offset = np.linalg.norm(of_v)
                
   
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
        print ok
        ok=unsout.setArrayF(component,"mass",mass_out)
        print ok
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
