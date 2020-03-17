import numpy as np
<<<<<<< HEAD
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
    def __init__(self, file_path,getcen=True,**kwargs):
        # get them vars !!!!ONLY RAMSES FOR NOW
        self.file = file_path
        self._center_history = np.array([0.,0.,0.])
        self.p = nbe.Info_sniffer(file_path)
        virial = kwargs.get('virial',False)
        comov = kwargs.get('comov',False)
        clumps = kwargs.get('clumps',False)
        self._dms, self._sts, self._gss  = False, False, False
        self.cen_done = False
        halo_vel = kwargs.get('halo_vel',[0.,0.,0.])    ##########
        # Load Data from Simulation
        self.n_tot, self.n_dm, self.n_st = nbe.check_particles(file_path)
        if self.n_dm > 0:
            print "loading Dark matter.."
            self.dm = _dark_matter(file_path, self.p,comov=comov)
            self._dms = True
        if self.n_st > 0:
            print "loading Stars.."
            self.st = _stars(file_path, self.p, comov=comov)
            self._sts = True
            print "loading Gas.."
            self.gs = _gas(file_path, self.p,comov=comov,virial=virial)
            self._gss = True
        else:    
                # if only DM is loaded computes center of zoom region 
                zoom_reg = self.dm.mass==self.dm.mass.min()
                dmcenter = nbe.real_center(self.dm.pos3d[zoom_reg],self.dm.mass[zoom_reg])
                self.center_shift(dmcenter)
                self.cen_done = True
        if (getcen) and not (self.cen_done):
            try:
                if (clumps):
                    # computes center with higher resolved clump
                    cen = self.dm.Clumps.pos3d[self.dm.Clumps.cell==self.dm.Clumps.cell.max()]
                else:
                    cen = nbe.real_center(self.st.pos3d,self.st.mass)
                self.center_shift(cen)
                self.cen_done = True
             
            except:
                print "No center cuz No clumps.. \nDo it yourself"    
           
    def r_virial(self, r_max,r_min=0,rotate=True,n=2.5, bins=5012):
        """
        The Function: main function of this story, calculates the three more commun
        virial radii: R200,R97 and R500. Also if rotate=True then it calculates the
        principal axes of the de stars structures and alines the hole box acordingly
        r_max: float, the maximal radius to prove in the search of the virial radii
        r_min: float, the minimal radius to prove in the search of the virial radii
        rotate: bool, if true rotates stars structure
        n: float, all particles with r>n*r200 will be ignore from now own.
        bins: int, number of bins for the coarse calculation
        """
=======
from . import  _dark_matter as d
from . import _stars as s
from . import _gas as g
import glob
from . import nbody_essentials as nbe


################################################################################
class Galaxy_Hound:
    def __init__(self, file_path,getcen=False,**kwargs):
        """
        Main Object: creates and object that loads all galaxy components and
        gives axes for numpy arrays for every variable from every particle 
        or cell.
        """
        # save path to snapshot
        self.file = file_path
        halonu = len(glob.glob(file_path+"/halo*"))
        if halonu==0:newage = False
        else:newage = True
        # initialize simu parameters
        self.p = nbe.Info_sniffer(file_path, newage=newage)
        # dicern on ramses versions

        ############################### ARGUMENTS ##############################
        gas             = kwargs.get('gas'     ,True ) # Load gas
        self.rmax_rot   = kwargs.get('rmax_rot',10   ) # rmax for rotation
        self.quiet      = kwargs.get('quiet'   ,False) # avoid printing
        hsml            = kwargs.get('hsml'    ,False) # force a res level
        self.dmo        = kwargs.get('dmo'     ,False) # dark matter only
        self.flush      = kwargs.get('flush'   ,False) # !!!!!!!!!!!!!!!! Check
        dens            = kwargs.get('dens'    ,False) # compute densities
        comov           = kwargs.get('comov'   ,True ) # comoving coordinates
        ########################################################################
        # flags for components
        self._dms, self._sts, self._gss  = False, False, False
        # inizialize mean halo velocity
        halo_vel = kwargs.get('halo_vel',[0.,0.,0.])    ########
        ############################### load data ##############################
        # read headers
        self.n_tot, self.n_dm, self.n_st = nbe.check_particles(file_path)
        # dark matter
        if self.n_dm > 0:
            if not self.quiet: print("loading Dark matter..")
            self.dm = d._dark_matter(file_path,self.p,comov=comov)
            self._dms = True
        # stars
        if self.n_st > 0 and self.dmo==False:
            if not self.quiet: print("loading Stars..")
            self.st = s._stars(file_path,self.p, comov=comov)
            self._sts = True
            # where there is stars there is gas
            if  gas==True:
                if not self.quiet: print("loading Gas..")
                self.gs = g._gas(file_path, self.p,comov=comov)
                self._gss = True
        else:
            self.dmo = True

    def r_virial(self,r_max=600,r_min=0,rotate=True,n=2.5,bins=512):
        """
        once the center have been defined this function computes
        several virial radii for the galaxy:
        r_200 : the radius where the spherical density is 200 times
        the critical density of the universe
        r_97 : the radius where the spherical density is 97 times
        the critical density of the universe
        r_BN : virial radius as defined by Brian & Norman 1998
        """
        # initializing
        rmax_rot = self.rmax_rot
>>>>>>> python3Ver
        positions = np.array([], dtype=np.int64).reshape(0,3)
        masses = np.array([], dtype=np.int64)
        # stack available masses and positions 
        if (self._dms):
            positions = np.vstack([positions,self.dm.pos3d])
            masses = np.append(masses,self.dm.mass)
        if (self._sts):
            positions = np.vstack([positions,self.st.pos3d])
            masses = np.append(masses,self.st.mass)
        if (self._gss):
            positions = np.vstack([positions,self.gs.pos3d])
            masses = np.append(masses,self.gs.mass)
<<<<<<< HEAD

        r = np.sqrt((positions[:,0])**2 +(positions[:,1])**2 +(positions[:,2])**2 )
        
        try:
            # calculating r200
            mhist, rhist = np.histogram(r,range=(0.0,r_max),bins=bins, weights=masses )
            vol_bin = (4./3.)*np.pi*(rhist[:-1]**3)
            r_bin = rhist[:-1]+ 0.5*(rhist[2]-rhist[1])
            rho_s = np.cumsum(mhist) / vol_bin
            # coarse calculation
            bin200 = np.argmin(np.abs(rho_s - (200 * self.p.rho_crit)))
            bin97 = np.argmin(np.abs(rho_s - (97 * self.p.rho_crit)))
            bin500 = np.argmin(np.abs(rho_s - (500 * self.p.rho_crit)))
            # fine calculation R200
            m_in = np.sum(masses[np.where(r<r_bin[bin200-1])])
            msum = np.cumsum(masses[np.where((r>=r_bin[bin200-1])&(r<r_bin[bin200+1]))])
            msum += m_in
            r_slice = r[np.where((r>=r_bin[bin200-1])&(r<r_bin[bin200+1]))] 
            vol_slice =  (4./3.)*np.pi*(r_slice**3)
            rho_slice =  msum / vol_slice
            self.r200= r_slice[np.argmin(np.abs(rho_slice-(200.*self.p.rho_crit)))]
            # fine calculation R97
            m_in = np.sum(masses[np.where(r<r_bin[bin97-1])])
            msum = np.cumsum(masses[np.where((r>=r_bin[bin97-1])&(r<r_bin[bin97+1]))])
            msum += m_in
            r_slice = r[np.where((r>=r_bin[bin97-1])&(r<r_bin[bin97+1]))] 
            vol_slice =  (4./3.)*np.pi*(r_slice**3)
            rho_slice =  msum / vol_slice
            self.r97= r_slice[np.argmin(np.abs(rho_slice-(97.*self.p.rho_crit)))]
            # fine calculation R500
            m_in = np.sum(masses[np.where(r<r_bin[bin500-1])])
            msum = np.cumsum(masses[np.where((r>=r_bin[bin500-1])&(r<r_bin[bin500+1]))])
            msum += m_in
            r_slice = r[np.where((r>=r_bin[bin500-1])&(r<r_bin[bin500+1]))] 
            vol_slice =  (4./3.)*np.pi*(r_slice**3)
            rho_slice =  msum / vol_slice
            self.r500= r_slice[np.argmin(np.abs(rho_slice-(500.*self.p.rho_crit)))]
            # marker 
            rnot = False
 
        except:
            sys.exit( "virial radius did not converged ")
        
        if (rotate)and((self._sts)or(self._gss)):
            print '| r_200 = {0:.3f} kpc'.format(self.r200)
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
    
=======
        # find virial radii 
        r = np.sqrt((positions[:,0])**2 +(positions[:,1])**2 +(positions[:,2])**2 )
        a,b,c,d = nbe.get_radii(r,masses,self.p,r_max,bins=bins)
        self.delta_crit, self.r200,self.r97,self.rBN = a,b,c,d
        rnot = False
        
        if (rotate)and(self._sts):
            if (self.flush):self.redefine(n,simple=True)
            if self.p.Z>2:
                self.rotate_galaxy(rmin=0.5,rmax=rmax_rot)
            else:
                self.rotate_galaxy()
            D = np.dot(self.matrix_T,np.dot(self.matrix_P,np.transpose(self.matrix_T)))
            if not self.quiet: nbe.print_matrix(D)

        if (rotate)and(self.dmo):
            if self.p.Z>2:
                self.rotate_galaxy(rmin=0.5,rmax=rmax_rot)
            else:
                self.rotate_galaxy(rmin=3,rmax=20,comp='dm')
            D = np.dot(self.matrix_T,np.dot(self.matrix_P,np.transpose(self.matrix_T)))
            if not self.quiet: nbe.print_matrix(D)

        self.frame_of_ref(r)
        self.redefine(n)

>>>>>>> python3Ver
    def center_shift(self,nucenter):
        """
        shifts the center of the whole box
        nucenter: [x,y,z] cordinates of the new center
        """
        self.center = np.zeros(3)
        self._center_history = np.vstack((self._center_history,nucenter))
        if (self._dms):
            self.dm.shift(nucenter)
        if (self._sts):
            self.st.shift(nucenter)
        if (self._gss):
            self.gs.shift(nucenter)
<<<<<<< HEAD
      
    def redefine(self,n):
        """
        cuts inside n*r200
        n: times r200 will be outer limit of the data after
        """
        if (self._dms):
            self.dm.halo_Only(self.center, n, self.r200)
        if (self._sts):
            self.st.halo_Only(self.center, n, self.r200, self.r97)
        if (self._gss):
            self.gs.halo_Only(self.center, n, self.r200)

    def rotate_galaxy(self,rmin=3,rmax=10):
        """
        rotates whole box according to the principal
        axes of the baryonic structure
        """
        if (self._sts):
            pos_ring = self.st.pos3d[(self.st.r<rmax)&(self.st.r>rmin)]
        else:
            pos_ring = self.gs.pos3d[(self.gs.r<rmax)&(self.gs.r>rmin)]
=======

    def redefine(self,n,simple=False):
        if (self._dms):
            self.dm.halo_Only(self.center, n, self.rBN,simple=simple)
        if (self._sts):
            self.st.halo_Only(self.center, n, self.rBN,simple=simple)
        if (self._gss):
            self.gs.halo_Only(self.center, n, self.rBN,simple=simple)

    def rotate_galaxy(self,rmin=3,rmax=10,comp='st',affect=True):
        if comp == 'st':
            r2 = (self.st.pos3d[:,0])**2 +(self.st.pos3d[:,1])**2 +(self.st.pos3d[:,2])**2
            pos_ring = self.st.pos3d[(r2<rmax**2)&(r2>rmin**2)]
        elif comp== 'dm':
            r2 = (self.dm.pos3d[:,0])**2 +(self.dm.pos3d[:,1])**2 +(self.dm.pos3d[:,2])**2
            pos_ring = self.dm.pos3d[(r2<rmax**2)&(r2>rmin**2)]

>>>>>>> python3Ver
        P = np.zeros((3,3))
        for i in range(3):
            for j in range(3):
                first = np.mean(pos_ring[:,i]*pos_ring[:,j])
                second =(np.mean(pos_ring[:,i])*np.mean(pos_ring[:,j]))
                P[i][j] = first - second
        eigen_values,evecs = np.linalg.eig(P)
        order = np.argsort(abs(eigen_values))
        T = np.zeros((3,3))
        T[0],T[1],T[2] = evecs[:,order[2]],evecs[:,order[1]],evecs[:,order[0]]
        self.matrix_T = T
        self.matrix_P = P
        if (self._dms) and (affect):
            self.dm.rotate(T)
        if (self._sts) and (affect):
            self.st.rotate(T)
        if (self._gss) and (affect):
            self.gs.rotate(T)
<<<<<<< HEAD
                
   
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
   
=======

    def frame_of_ref(self,r):
        ############################################################
        # substract the average speed of the system
        positions = vels = np.array([], dtype=np.int64).reshape(0,3)
        mass = np.array([])
        if (self._dms):
            #r = np.append(r,self.dm.r)
            mass = np.append(mass,self.dm.mass)
            vels = np.vstack([vels,self.dm.vel3d])

        if (self._sts):
            #r = np.append(r,self.st.r)
            mass = np.append(mass,self.st.mass)
            vels = np.vstack([vels,self.st.vel3d])
        if (self._gss):
            #r = np.append(r,self.gs.r)
            mass = np.append(mass,self.gs.mass)
            vels = np.vstack([vels,self.gs.vel3d])
        # velocity of the center of mass 
        sel = np.where(r<self.r200)
        self.com_vx = np.sum(mass[sel]*vels[sel,0])/np.sum(mass[sel])
        self.com_vy = np.sum(mass[sel]*vels[sel,1])/np.sum(mass[sel])
        self.com_vz = np.sum(mass[sel]*vels[sel,2])/np.sum(mass[sel])
        if (self._dms):self.dm.vel_frame(self.com_vx,self.com_vy,self.com_vz)
        if (self._sts):self.st.vel_frame(self.com_vx,self.com_vy,self.com_vz)
        if (self._gss):self.gs.vel_frame(self.com_vx,self.com_vy,self.com_vz)
        ##########################################################


>>>>>>> python3Ver
