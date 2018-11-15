import numpy as np


class HALOBHydro:
    def __init__(self,where="home", **kargs):
        self.name = "Halo B"
        self.namenospace = "HALOB"
        self.version = "hydro" 
        if where=="manu":
            self.path = "/net/direct/backup/pol/Box20Mpc-Zoom-MWGalax-MUSIC-bis/Halo24/Zoom4-Halo24-refmap128-HydroRun/output_00417"
        else:    
            self.path = "/data/POL/HALOB/hydro/output_00417"
        self.c_dm_com = np.array([ 9667.68750, 9866.01957, 9801.56035])# kpc
        self.c_rho_st = np.array([ 9666.932  , 9864.464  , 9799.166])# kpc
        self.rmax = 794.76 #kpc
        self.r200 = 177.54 #kpc
        self.dm_part_mass = 2.30812e5 #Msun
        self.pot_max = -29559.54 #pc km^2 M_sun s^-2
        self.potfile = "/home/anunez/WKBL_candidates/Catalog/potentials/Psi_halo_B_DM_baryons_Rmax=805.0kpc_dimensionful.txt"
        self.q = -0.11878
        self.M200_dm =  5.00696645632e+11 #Msun
        self.M200_st = 0.79612198912e+11 #Msun
        self.M10_st   = 7.50630e+10 #Msun
        self.Mfire_st = 4.93454e+10 #Msun
        self.dm_fit  = np.array([7.622, 4.885, 2.639, 2.621, 0.127])
        self.st_fit   = np.array([8.095,  0.340, 0.260, 2.297, 0.540, 8.634, 0.608, 3.104, 7.885, 1.219, 6.159])
        self.gs_fit   = np.array([9.037 , 0.031, 6.116])
        #v_esc for the 4 testing points vs Eddington
        # spherical
        self.v_esc_sph       = np.array([567.60, 497.32, 384.80, 303.38, 198.36 ])
        self.v_esc_sigma_sph = np.array([12.57, 10.09, 5.89, 1.93, 12.19])
        # Isopotential
        self.v_esc_iso       = np.array([596.83, 510.02, 395.56, 291.51, 191.31 ])
        self.v_esc_sigma_iso = np.array([9.48, 8.17, 5.33, 2.71, 3.79 ])




class HALOBdmo:
    def __init__(self,where="home", **kargs):
        self.name = "Halo B"
        self.namenospace = "HALOB"
        self.version = "DMO" 
        if where=="manu":
            self.path = "/net/direct/backup/pol/Box20Mpc-Zoom-MWGalax-MUSIC-bis/Halo24/Zoom4-Halo24-refmap128-HydroRun/output_00417"
        else:    
            self.path = '/data/POL/HALOB/DMO/output_00041'
        self.c_dm_com = np.array([ 9693.93650, 9871.29969, 9805.710 ])# kpc
        self.c_dm_rho = np.array([9692.29980469, 9869.86230469, 9804.25097656])# kpc
        self.c_dm_pot = np.array([9692.24414062, 9869.58789062, 9804.03222656])# kpc
        self.rmax = 787.79 #kpc
        self.r200 = 163.28 #kpc
        self.dm_part_mass = 2.75776e5 #Msun
        self.pot_max = -28066.64 #pc km^2 M_sun s^-2
        self.potfile = "/home/anunez/WKBL_candidates/Catalog/potentials/Psi_halo_B_DMO_Rmax=805.0kpc_dimensionful.txt" 
        self.q = -0.18997 
        self.M200_dm =  4.960e+11 #Msun
        self.dm_fit  = [7.101,10.338,1,2.785,1.101]
        #v_esc for the 4 testing points vs Eddington
        # spherical
        self.v_esc_sph       = np.array([368.25, 352.71, 318.96, 261.69, 196.46  ])
        self.v_esc_sigma_sph = np.array([10.61, 7.29, 6.63, 3.86, 5.11 ])
        # Isopotential
        self.v_esc_iso       = np.array([400.91, 371.06, 323.31, 259.02, 171.1])
        self.v_esc_sigma_iso = np.array([7.72, 4.97, 5.32, 4.69, 5.32 ])



class HALOCHydro:
    def __init__(self,where="home", **kargs):
        self.namenospace = "HALOC"
        self.name = "Halo C"
        self.version = "hydro"
        if where=="manu":
            self.path = "/net/direct/backup/pol/Box20Mpc-Zoom-MWGalax-MUSIC-bis/Halo19/Zoom4-Halo19-refmap128-Hydro/output_00440"
        else:
            self.path = "/data/POL/HALOC_19/Hydro/output_00440"
        self.pot_max  = -21926.61  #pc km^2 M_sun s^-2
        self.rmax = 1725.09 #kpc
        self.c_dm_com = np.array([ 9868.14825,  9745.19211,  9766.90090 ])# kpc
        self.c_rho_st = np.array([ 9867.712,    9743.975,    9765.676 ])# kpc
        self.q        = 0.057
        self.r200     = 182.23 #kpc
        self.M200_dm  = 5.50e+11 #Msun
        self.M200_st = 9.56e+10 #Msun
        self.M10_st   = 8.91017e+10 #Msun
        self.Mfire_st = 2.99687e+10 #Msun
        self.dm_fit  = np.array([7.746, 4.210, 2.165, 2.504, 0.189])
        self.st_fit   = np.array([9.581, 0.129, 0.818, 4.896, 0.330, 9.080, 0.021, 1.622, 7.860, 0.674, 8.671])
        self.gs_fit   = np.array([9.327, 0.135, 1.223])
        #v_esc for the 4 testing points vs Eddington
        # spherical
        self.v_esc_sph       = np.array([601.24, 514.17, 403.19, 332.32, 234.95 ])
        self.v_esc_sigma_sph = np.array([14.03, 7.52, 7.89, 8.91, 16.47 ])
        # Isopotential
        self.v_esc_iso       = np.array([635.83, 531.29, 418.26, 323.69, 236.91 ])
        self.v_esc_sigma_iso = np.array([9.00, 6.28, 3.67, 2.46, 11.43 ])



class HALOCdmo:
    def __init__(self,where="home", **kargs):
        self.name = "Halo C"
        self.namenospace = "HALOC"
        self.version = "DMO"
        if where=="manu":
            self.path = "/net/direct/backup/pol/Box20Mpc-Zoom-MWGalax-MUSIC-bis/Halo19/Zoom4-Halo19-refmap128-DMonly/output_00041"
        else:
            self.path = "/data/POL/HALOC_19/DMO/output_00041"
        self.pot_max  = -20588.90  #pc km^2 M_sun s^-2
        self.rmax = 1688.46 #kpc
        self.c_dm_com = np.array([ 9861.9777, 9776.37011, 9790.41548 ])# kpc
        self.c_dm_pot= np.array([9860.93359375, 9775.50292969, 9789.64746094])# kpc
        self.q = -0.19327
        self.M200_dm =  6.248e+11 #Msun
        self.dm_fit  = rho = [6.848,14.291,1,2.805,1.065]
        #v_esc for the 4 testing points vs Eddington
        # spherical
        self.v_esc_sph       = np.array([368.98, 367.24, 337.66, 283.03, 218.01 ])
        self.v_esc_sigma_sph = np.array([11.16, 5.56, 3.88, 6.15, 11.10 ])
        # Isopotential
        self.v_esc_iso       = np.array([411.24, 394.22, 354.46, 328.61, 276.17 ])
        self.v_esc_sigma_iso = np.array([7.54, 4.38, 4.84, 6.55, 51.55])






class MochimaHydro:
    def __init__(self,where="home", **kargs):
        self.name = "Mochima"
        self.namenospace = self.name 
        self.version = "hydro"
        if where=="manu":
            self.path = " "
        else:
            self.path = "/data/OWN/SF1test/SF0/mstar1_T3600/output_00041"
        self.pot_max  = -27342.01  #pc km^2 M_sun s^-2
        self.rmax = 2908.43 #kpc
        self.c_dm_com = np.array([20418.88714, 17567.72033, 17124.40448 ])# kpc
        self.c_rho_st = np.array([20415.127, 17564.615, 17121.342])# kpc
        self.q = -0.003443
        self.M200_dm  = 8.21482749952e+11 #Msun
        self.M200_st  = 1.17450924032e+11 #Msun
        self.M10_st   = 1.08706e+11 #Msun
        self.Mfire_st = 2.14903e+08 #Msun
        self.dm_fit  = np.array([7.414, 7.364, 1.768, 2.633, 0.468])
        self.st_fit   = np.array([9.381, 0.678, 0.381, 1.063, 0.677, 7.384, 1.049, 2.098, 8.892, 1.067, 3.017])
        self.gs_fit   = np.array([8.497, 0.119, 5.971])
        #v_esc for the 4 testing points vs Eddington
        # spherical
        self.v_esc_sph       = np.array([612.71, 552.56, 446.87, 371.59, 259.18 ])
        self.v_esc_sigma_sph = np.array([24.56, 15.05, 8.67, 9.39, 5.93 ])
        # Isopotential
        self.v_esc_iso       = np.array([693.02, 594.51, 480.69, 371.71, 259.58 ])
        self.v_esc_sigma_iso = np.array([9.40, 7.34, 4.26, 2.70, 2.14 ])




class Mochimadmo:
    def __init__(self,where="home", **kargs):
        self.name = "Mochima"
        self.namenospace = self.name 
        self.version = "DMO"
        if where=="manu":
            self.path = " "
        else:
            self.path = "/data/OWN/DMO/mochima2_Z5/output_00041"
        self.pot_max  = -28806.34  #pc km^2 M_sun s^-2
        self.rmax = 2792.98 #kpc
        self.c_dm_com = np.array([ 20438.06143, 17580.66710, 17120.75325])# kpc
        self.c_dm_pot = np.array([20432.74609, 17575.785525, 17116.07815])# kpc
        self.q = -0.16956
        self.M200_dm =  9.13551392768e+11 #Msun
        self.dm_fit  = rho = [6.963,13.786,1,2.721,1.066]
        #v_esc for the 4 testing points vs Eddington
        # spherical
        self.v_esc_sph       = np.array([381.01, 396.13, 378.41, 333.66, 270.35 ])
        self.v_esc_sigma_sph = np.array([31.45, 20.37, 9.69, 7.92, 5.47])
        # Isopotential
        self.v_esc_iso       = np.array([473.15, 455.04, 410.66, 347.44, 256.42 ])
        self.v_esc_sigma_iso = np.array([8.24, 4.91, 4.49, 4.40, 3.33 ])





class AdicoraHydro:
    def __init__(self,where="home", **kargs):
        self.name = "Adicora"
        self.namenospace = self.name 
        self.version = "hydro"
        if where=="manu":
            self.path = " "
        else:
            self.path = "/data/OWN/Adicora/SF0/Stable/output_00041"
        self.pot_max  = -24172.00 #pc km^2 M_sun s^-2
        self.rmax = 1580# rmax +130, rmax_true = 1453.96 #kpc
        self.c_dm_com = np.array([14313.79353149, 15227.11111001, 15695.32341304 ])# kpc
        self.c_rho_st = np.array([14314.327, 15226.841, 15695.589])# kpc
        self.q = -0.19522
        self.M200_dm  = 8.10753261568e+11 #Msun
        self.M200_st  = 1.10771314688e+11 #Msun
        self.M10_st   = 1.01098e+11 #Msun
        self.Mfire_st = 6.36326e+09 #Msun
        self.dm_fit  = np.array([7.6676,6.6669,1.4075,2.7560,0.0898])
        self.st_fit   = np.array([9.5591,0.6407,3.7140,0.6177,0.9843,9.3026,0.7869,0.8915,8.4928,0.9000,4.8636])
        self.gs_fit   = np.array([8.415 , 0.131,7.510])
        #v_esc for the 4 testing points vs Eddington
        # spherical
        self.v_esc_sph       = np.array([612.71, 552.56, 446.87, 371.59, 259.18 ])
        self.v_esc_sigma_sph = np.array([24.56, 15.05, 8.67, 9.39, 5.93])
        # Isopotential
        self.v_esc_iso       = np.array([691.22, 592.41, 478.10, 368.35, 254.75  ])
        self.v_esc_sigma_iso = np.array([9.43, 7.37, 4.28, 2.72, 2.18  ])



class Adicoradmo:
    def __init__(self,where="home", **kargs):
        self.name = "Adicora"
        self.namenospace = self.name 
        self.version = "DMO"
        if where=="manu":
            self.path = " "
        else:
            self.path = "/data/OWN/DMO/Adicora/output_00041"
        self.pot_max  = -23323.76  #pc km^2 M_sun s^-2
        self.rmax =  1580# rmax + 130, rmax_true = 1447.53 #kpc
        self.c_dm_com = np.array([14308.07419, 15223.16481, 15686.17544]) #kpc
        self.c_dm_pot = np.array([14307.2343, 15221.832031, 15684.78808]) #kpc
        self.q = -0.15649
        self.M200_dm = 9.16783366144e+11 #Msun
        self.dm_fit  = rho = [7.1813,11.2881,1.0000,2.7440,0.9748]
        #v_esc for the 4 testing points vs Eddington
        # spherical
        self.v_esc_sph       = np.array([427.74, 426.13, 391.99, 329.54, 250.96 ])
        self.v_esc_sigma_sph = np.array([14.88, 9.06, 6.05, 3.95, 3.76])
        # Isopotential
        self.v_esc_iso       = np.array([475.25, 452.81, 404.40, 337.74, 235.6  ])
        self.v_esc_sigma_iso = np.array([8.33, 5.23, 3.91, 3.12, 3.69])




