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
        self.unitl=_vars["unit_l"]/(3.08e18)*self.pctocm
        self.unitd=_vars["unit_d"]/(self.pctocm/(3.08*10**18))**3
        self.unitt=_vars["unit_t"]
        self.simutokpc=self.unitl*_vars["H0"]/self.pctocm/1e5
        self.simutoMsun=(self.unitd*self.unitl**3)/1000/self.msuntokg
        self.unitsimutoMsunkpc3=self.unitd*self.pctocm**3/1000/self.msuntokg
        self.simutokms = self.unitl/10**5/self.unitt
        self.kgtoGeV = 1/1.783e-27
        self.kpctocm = 3.086e21
        self.G = 6.67384e-11 * self.msuntokg / (3.08567758e19**3)
        self.rho_crit = (3 * (self.H0**2) / 3.08567758e19**2)/ (8*np.pi*self.G)
        self.simutoGeVcm3 = (self.simutoMsun*self.msuntokg*self.kgtoGeV) / (self.simutokpc*self.kpctocm)**3

def _get_center(output,clumps=False):
        """
        gets center of more resolved halo
        taken from hast library witten by V.perret
        https://bitbucket.org/vperret/hast 
        """
        p = Info_sniffer(output)
        list = glob.glob(output+'/clump_?????.txt?????')
        i=0
        for file in list:
                data = np.loadtxt(file,skiprows=1,dtype=None)
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
        if not (clumps):
            return data_sorted[0,4:7]
        else:
            return data_all[array.argsort()]

def get_com(pos,m):
    return np.array([np.dot(pos[:,0],m),
                     np.dot(pos[:,1],m),
                     np.dot(pos[:,2],m)])/ np.sum(m)
def r(pos):
    return np.sqrt(pos[:,1]**2 + pos[:,1]**2 + pos[:,2]**2)


def real_center(pos, m):
    """
    this method computes the center of mass of a recursively reducing
    sphere center in the previous COM as descrived by Schaller et al. 
    in
    arxiv:1505.05470v2
    input: dm, st, gas object from Galaxy_hound()
    coordinates of the real center of mass
 
    """
    final = np.zeros((1,3))
    while len(pos) > 200:
        com = get_com(pos,m)
        final += com
        pos -= com
        m = m[(r(pos)<0.9*np.max(r(pos)))]
        pos = pos[(r(pos)<0.9*np.max(r(pos)))]
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
    m_tot, m_sort = mass.sum() / 2., mass[np.argsort(r)]
    aux,counter = 0,-1
    while aux < m_tot:
        counter += 1
        aux += m_sort[counter]
    return r[counter]

