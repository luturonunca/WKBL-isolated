import configparser

class dm_input():
    def parse_input(self, input_file):
        inputs = configparser.SafeConfigParser()
        inputs.read(input_file)
        try:
            open(input_file,'r')
            print "### using  ",input_file," as inputs"
        except:
            print "### no input file named: ", input_file,", using default values"
        try:
            self.v_Sun = inputs.getfloat('general','v_sun')
        except:
            self.v_Sun = 220.
        try:
            self.v_esc = inputs.getfloat('general','v_esc')
        except:
            self.v_esc = 544.
        try:
            self.v_esc2 = inputs.getfloat('general','v_esc2')
        except:
            self.v_esc2 = 394.
        try:
            self.v_esc3 = inputs.getfloat('general','v_esc3')
        except:
            self.v_esc3 = 294.
        try:
            self.v_shm = inputs.getfloat('general','v_0shm')
        except:
            self.v_shm = 270.
        try:
            self.v_lin = inputs.getfloat('general','v_0lin')
        except:
            self.v_lin = 544.
        try:
            self.v_mao = inputs.getfloat('general','v_0mao') * self.v_esc
        except:
            self.v_mao = 0.13 * self.v_esc
        try:
            self.rho = inputs.getfloat('general','rho')
        except:
            self.rho = 0.3
        try:
            self.sigma_sd = inputs.getfloat('general','sigma_sd') / 0.3894
        except:
            self.v_esc = 2.19e-6 / 0.3894
        try:
            self.sigma_si = inputs.getfloat('general','sigma_si') / 0.3894
        except:
            self.v_esc = 2.85e-9 / 0.3894
        try:
            self.m_p = inputs.getfloat('general','m_p')
        except:
            self.m_p = 938.272e-3
        self.r_h = ((0.9*(self.m_p**(1/3)))+ 0.3) * (10**(-13)) # nuclear radius
        self.n_i = 0.772 ####CHECK!!!### # nuclear number density in the Sun
        self.v_0 =270.1e5 # velocity dispertion 
        self.h = 6.582e-25 # plank constant in GeV*s

