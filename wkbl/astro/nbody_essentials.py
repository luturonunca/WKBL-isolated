import sys
import math
import glob
import cmath
#import emcee
import subprocess
import numpy as np
from py_unsio import *
import scipy.special as sp
from numpy import exp, sqrt
import scipy.integrate as integrate
from sklearn.neighbors import KDTree
import warnings
warnings.filterwarnings('ignore')


class Info_sniffer:
    def __init__(self, file_path):
        _vars=dict()
        file = glob.glob(file_path+"/info_?????.txt")[0]
        with open(file) as outinfo:
            for line in outinfo:
                eq_index = line.find('=')
                if eq_index == -1:
                    continue
                var_name = line[:eq_index].strip()
                number = float(line[eq_index + 1:].strip())
                _vars[var_name] = number
                if var_name == "unit_t":
                    break
        self.H0=_vars["H0"]
        self.aexp=_vars["aexp"]
        self.msuntokg = 1.99844e30
        self.pctocm = 3.08567758e18
        self.cmtopc = 1/self.pctocm
        
        self.unitl=_vars["unit_l"]#/(3.08e18)*self.pctocm
        self.unitd=_vars["unit_d"]/(self.pctocm/(3.08*10**18))**3
        self.unitt=_vars["unit_t"]
        self.simutokpc=self.unitl*_vars["H0"]/self.pctocm/1e5
        self.simutocm=self.unitl*_vars["H0"]/1e5
        self.simutoMsun=(self.unitd*self.unitl**3)/1000/self.msuntokg
        self.unitsimutoMsunkpc3=self.unitd*self.pctocm**3/1000/self.msuntokg
        self.simutokms = self.unitl/10**5/self.unitt
        self.kgtoGeV = 1/1.783e-27
        self.kpctocm = 3.086e21
        self.k_boltz = 1.3806488e-23 * 1e-6 / self.msuntokg # Msun * km**2 / s**2
        self.m_p = 1.67262177e-27/ self.msuntokg # Msun
        self.Z = (1./self.aexp) - 1.
        self.kpctokm = self.kpctocm / 1e5
        self.G = 6.67384e-11 * self.msuntokg / (3.08567758e19**3)
        self.rho_crit = (3 * (self.H0**2) / 3.08567758e19**2)/ (8*np.pi*self.G)
        self.simutoGeVcm3 = (self.simutoMsun*self.msuntokg*self.kgtoGeV) / (self.simutokpc*self.kpctocm)**3
        self.mH = 1.6600000e-24
        self.kB = 1.3806200e-16
        self.simutoKelvin =  self.mH/self.kB * (self.simutocm/self.unitt)**2

def _get_center(output,clumps=False, sf_hist=False):
        """
        gets center of more resolved halo
        taken from hast library witten by V.perret
        https://bitbucket.org/vperret/hast 
        """
        p = Info_sniffer(output)
        if (sf_hist):
            list = glob.glob(output+'/stars_?????.out?????')
        else:
            list = glob.glob(output+'/clump_?????.txt?????')
        i=0
        for file in list:
                try:
                    data = np.loadtxt(file,skiprows=1,dtype=None)
                except:
                    continue
                if(np.size(data)==0):
                        continue
                if(i>0):
                        data_all = np.vstack((data_all,data))
                else:
                        data_all = data
                i=i+1
        array=(1e4*data_all[:,3]/np.max(data_all[:,3]))*(data_all[:,8]/np.max(data_all[:,8])).astype(int, copy=False)
        data_sorted = data_all[array.argsort()]
        data_sorted = data_sorted[::-1]
        if (sf_hist):
            return data_all
        elif not (clumps):
            return data_sorted[0,4:7]
        else:
            return data_all[array.argsort()]

def get_com(pos,m):
    return np.array([np.dot(pos[:,0],m),
                     np.dot(pos[:,1],m),
                     np.dot(pos[:,2],m)])/ np.sum(m)
def get_r(pos):
    return np.sqrt(pos[:,1]**2 + pos[:,1]**2 + pos[:,2]**2)


def real_center(pos, mass, n=7000):
    """
    this method computes the center of mass of a recursively reducing
    sphere center in the previous COM as descrived by Schaller et al. 
    in
    arxiv:1505.05470v2
    input: dm, st, gas object from Galaxy_hound()
    coordinates of the real center of mass
 
    """
    p = np.copy(pos)
    m = np.copy(mass)
    final = np.zeros((1,3))
    while len(p) > n:
        com = get_com(p,m)
        final += com
        p -= com
        r = get_r(p)
        m = m[np.where(r<0.9*r.max())]
        p = p[np.where(r<0.9*r.max())]
    return final[0]

def matrix_vs_vector(mat,vec):
    """
    costume matrix vector multiplication to speed up the product between a matrix
    rotation matrix for example, and an array of vectors
    mat --> matrix
    vec --> array of 3d vectors
    output: the resulting array of 3d vectors
    """
    res = np.zeros((len(vec),3))
    res[:,0] = mat[0][0]*vec[:,0] + mat[0][1]*vec[:,1] + mat[0][2]*vec[:,2]
    res[:,1] = mat[1][0]*vec[:,0] + mat[1][1]*vec[:,1] + mat[1][2]*vec[:,2]
    res[:,2] = mat[2][0]*vec[:,0] + mat[2][1]*vec[:,1] + mat[2][2]*vec[:,2]
    return res


def all_inside(pos3d, center, r_search):
    try:
        tree = KDTree(np.squeeze(pos3d))
        in_halo = tree.query_radius(center,r_search)[0]
    except:
        sys.exit("Nope")
    return in_halo

def half_mass(mass,r):
    """
    returns half mass radius
    i.e the radius at which the half of the total 
    mass is contained
    input:
    mass = array of particles mass
    r =  array of particles distance from center 
    """
    m_tot, m_sort = mass.sum() / 2., mass[np.argsort(r)]
    aux,counter = 0,-1
    while aux < m_tot:
        counter += 1
        aux += m_sort[counter]
    return r[counter]


def FIRE_st_mass(mass,r,r97):
    """
    from arxiv:1702.06148v1 footnote 9
    returns the mass defined there as M_st
    and the radius they define
    to plot on top of the Moster band
    """
    R15 = 0.15*r97
    R_half = half_mass(mass[r< R15],r[r< R15])
    R_FIRE = half_mass(mass[r< 3*R_half],r[r< 3*R_half])
    return mass[r < R_FIRE].sum(), R_FIRE
    

def check_particles(path):
    """
    this routine check the particle number in  the simulation
    by reading the header file in the snapshot file and using
    as input the path to the stapshot.
    """
    list = glob.glob(path+"/header_?????.txt")
    linum = n = 0
    part_arrays = fi = np.array([])
    for f in list:
        fi = np.append(fi,f)
        n+=1
    myfi = open(fi[0])

    for l in myfi:
        linum += 1
        if linum%2==0 and linum<=8:
            row = l.split(' ')
            part_arrays = np.append(part_arrays,np.float(row[-1]))
    n_tot = part_arrays[0]
    n_dm = part_arrays[1]
    n_st = part_arrays[2]
    if (n_tot>0):
        return n_tot, n_dm, n_st
    else:
        myfi, linum = open(fi[0]), -1
        for l in myfi:
            linum += 1
            row = l.split('  ')
            if linum==7:
                n_dm = float(row[-1])
            if linum==8:
                n_st = float(row[-1])
        n_tot = n_st + n_dm
        return int(n_tot), int(n_dm), int(n_st)

