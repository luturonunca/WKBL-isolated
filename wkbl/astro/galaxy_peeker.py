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

#### master branch ####
 
 

class Galaxy_Hound:
    def __init__(self, file_path, component,getcen=True,**kwargs):
        # get them vars !!!!ONLY RAMSES FOR NOW
        self.file = file_path
        self._center_history = np.array([0.,0.,0.])
        self.p = nbe.Info_sniffer(file_path)
        comov = kwargs.get('comov',False)
        self._dms, self._sts, self._gss  = False, False, False
        self.cen_done = False
        halo_vel = kwargs.get('halo_vel',[0.,0.,0.])    ##########
        # Load Data from Simulation
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
            if component=="halo":
                # if only DM is loaded computes center of zoom region               
                zoom_reg = self.dm.mass==self.dm.mass.min()
                dmcenter = nbe.real_center(self.dm.pos3d[zoom_reg],self.dm.mass[zoom_reg])
                self.center_shift(dmcenter)
                self.cen_done = True
        if (self._sts):
            print "loading Stars.."
            self.st = _stars(file_path, self.p, comov=comov)
        if (self._gss):
            print "loading Gas.."
            self.gs = _gas(file_path, self.p,comov=comov)
        if (getcen) and not (self.cen_done):
            try:
                # computes center with higher resolved clump
                cen = self.dm.Clumps.pos3d[self.dm.Clumps.cell==self.dm.Clumps.cell.max()]
                self.center_shift(cen)
            except:
                print "No center cuz No clumps.. \nDo it yourself"    
           
    def r_virial(self, r_max,r_min=0,rotate=True,n=2.5, bins=512):
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
        
        #try:
        for i in range(1):
            mhist, rhist = np.histogram(r,range=(0.0,r_max),bins=bins, weights=masses )
            vol_bin = (4./3.)*np.pi*(rhist[:-1]**3)
            r_bin = rhist[:-1]+ 0.5*(rhist[2]-rhist[1])
            rho_s = np.cumsum(mhist) / vol_bin
            self.r200 = r_bin[np.argmin(np.abs(rho_s - (200 * self.p.rho_crit)))]
            self.r97 = r_bin[np.argmin(np.abs(rho_s - (97 * self.p.rho_crit)))]
            rnot = False
 
        #except:
        #    print "virial radius did not converged "
        #    sys.exit()
        if (rotate)and((self._sts)or(self._gss)):
            print '| r_200 = {0}'.format(self.r200)
            print '---- taking particles inside {0} * r200'.format(n)
            self.redefine(n)
            print '| number of praticles inside {0} * r200 '.format(n)
            if (self._dms):
                print '| dm mass       =  {0:1.3e} M_sun'.format(self.dm.total_m)
                print '| p_dm_200      =  {0:1.3e} particles'.format(len(self.dm.pos3d))
            if (self._sts):
                print '| stellar mass  =  {0:1.3e} M_sun'.format(self.st.total_m)
                print '| p_st_200      =  {0:1.3e} psrticles'.format(len(self.st.pos3d))
            if (self._gss):
                print '| gas mass      =  {0:1.3e} M_sun'.format(self.gs.total_m)
                print '| p_gs_200      =  {0:1.3e} particles'.format(len(self.gs.pos3d))
            print '---- rotating galaxy '
            self.rotate_galaxy()
            D = np.dot(self.matrix_T,np.dot(self.matrix_P,np.transpose(self.matrix_T)))
            print '| Diagonal matrix computed '
            print '|    |{0:2d},{1:2d},{2:2d}|'.format(int(D[0,0]),int(D[0,1]),int(D[0,2]))
            print '| D =|{0:2d},{1:2d},{2:2d}|'.format(int(D[1,0]),int(D[1,1]),int(D[1,2]))
            print '|    |{0:2d},{1:2d},{2:2d}|'.format(int(D[2,0]),int(D[2,1]),int(D[2,2]))
        elif (rotate):
            self.redefine(n)
    
    def center_shift(self,nucenter):
        self.center = np.zeros(3)
        self._center_history = np.vstack((self._center_history,nucenter))
        if (self._dms):
            self.dm.shift(nucenter)
        if (self._sts):
            self.st.shift(nucenter)
        if (self._gss):
            self.gs.shift(nucenter)
      
    def redefine(self,n):
        if (self._dms):
            self.dm.halo_Only(self.center, n, self.r200)
        if (self._sts):
            self.st.halo_Only(self.center, n, self.r200, self.r97)
        if (self._gss):
            self.gs.halo_Only(self.center, n, self.r200)

    def rotate_galaxy(self,rmin=3,rmax=10):
        if (self._sts):
            pos_ring = self.st.pos3d[(self.st.r<rmax)&(self.st.r>rmin)]
        else:
            pos_ring = self.gs.pos3d[(self.gs.r<rmax)&(self.gs.r>rmin)]
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
                
   
    def save_galaxy(self, name, fltype,density=False):
        unsout=CunsOut(name,fltype)
        ages = False
        if (self._dms):
            length = len(self.dm.pos3d)
            pos_out = self.dm.pos3d.reshape(length*3).astype(np.float32, copy=False)
            mass_out = self.dm.mass.astype(np.float32, copy=False)
            unsout.setArrayF("halo","mass",mass_out) # save mass
            unsout.setArrayF("halo","pos",pos_out)
            if (density):
                ok,rho_dm,_= CF.getDensity(np.array(self.dm.pos3d.reshape(len(self.dm.pos3d)*3),
                                           dtype=np.float32), self.dm.mass)
                if (ok):
                    rho_out_dm = rho_dm.astype(np.float32, copy=False)
                    unsout.setArrayF("halo","rho",rho_out_dm)
                    
        if (self._sts):
            length = len(self.st.pos3d)
            pos_out = self.st.pos3d.reshape(length*3).astype(np.float32, copy=False)
            mass_out = self.st.mass.astype(np.float32, copy=False)
            age_out = self.st.age.astype(np.float32, copy=False)
            ages = True
            unsout.setArrayF("stars","pos",pos_out)
            unsout.setArrayF("stars","mass",mass_out)
            unsout.setArrayF("stars","age",mass_out)
            if (density):
                ok,rho_st,_= CF.getDensity(np.array(self.st.pos3d.reshape(len(self.st.pos3d)*3),
                                           dtype=np.float32), self.st.mass)
                if (ok):
                    rho_out_st = rho_st.astype(np.float32, copy=False)
                    unsout.setArrayF("stars","rho",rho_out_st)
        if (self._gss):
            length = len(self.gs.pos3d)
            pos_out = self.gs.pos3d.reshape(length*3).astype(np.float32, copy=False)
            mass_out = self.gs.mass.astype(np.float32, copy=False)
            rho_out = self.gs.rho.astype(np.float32, copy=False)
            unsout.setArrayF("gas","pos",pos_out)
            unsout.setArrayF("gas","mass",mass_out)
            unsout.setArrayF("gas","rho",rho_out)
        unsout.save()
   
