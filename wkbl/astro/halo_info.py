import numpy as np


class HALOBHydro:
    def __init__(self,where="home", **kargs):
        self.name = "Halo B"
        self.version = "hydro" 
        if where=="manu":
            self.path = "/net/direct/backup/pol/Box20Mpc-Zoom-MWGalax-MUSIC-bis/Halo24/Zoom4-Halo24-refmap128-HydroRun/output_00417"
        else:    
            self.path = "/data/POL/HALOB/hydro/output_00417"
        self.c_st_rho = np.array([-0.06085646, -1.06476205, -1.4065722 ])# kpc
        self.c_st_com = np.array([ 0.0005449 ,  0.00235399, -0.00282367])# kpc
        self.c_dm_com = np.array([ 0.2661515 , -0.07751549,  0.41966864])# kpc
        self.rmax = 794.76 #kpc
        self.pot_max = -29616.08 #pc km^2 M_sun s^-2
        self.potfile = "/home/anunez/WKBL_candidates/Catalog/potentials/Psi_halo_B_DM_baryons_Rmax=805.0kpc_dimensionful.txt" 

class HALOBdmo:
    def __init__(self, **kargs):
        self.name = "Halo B"
        self.version = "DMO" 
        self.path = "/net/direct/backup/pol/Box20Mpc-Zoom-MWGalax-MUSIC-bis/Halo24/Zoom4-Halo24-refmap128-DMonly/output_00041"
        self.c_dm_com = np.array([ 0.03749188, 0.2631074,  0.2292752 ])# kpc
        self.rmax = 787.79 #kpc
        self.pot_max = -28064.58 #pc km^2 M_sun s^-2
        self.potfile = "/home/anunez/WKBL_candidates/Catalog/potentials/Psi_halo_B_DMO_Rmax=805.0kpc_dimensionful.txt" 

