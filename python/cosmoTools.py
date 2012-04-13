from scipy.integrate import quad
import numpy
        

class cosmology:
    """useful functions dependent upon cosmological parameters"""

    
    def __init__(self, omega_m0 = 0.274, omega_lam0 = 0.726, omega_r0 = None, h=0.71):
       
        self.omega_m0 = omega_m0
        self.omega_lam0 = omega_lam0
        if omega_r0 is None:
            self.omega_r0 = 4.15e-5/(h*h) # assumes 3 massless neutrino species
        else:
            self.omega_r0 = omega_r0
            
        self.h = h
        self.c = 2.99792458e5    # in km/s
        
        self.H0 = 100.0*h
        self.omega0 = self.omega_m0 + self.omega_lam0 + self.omega_r0
        
        # determine kappa = H0^2/c^2 * (omega0 - 1)
        self.kappa = self.H0**2/self.c**2 * (self.omega0 - 1.0)

    def E(self, z):
        """computes E(z)"""
        
        e = (self.omega_m0 *(1+z)**3 + self.omega_r0*(1+z)**4 + self.omega_lam0 + (1+self.omega0) )**(.5)

        return e

    def Hubble(self, z):
        """computes H(z)"""

        H_0 = 100.*self.h
        e = self.E(z)

        return H_0*e
    
    
    def comovingDist(self, z):
        """computes comoving distances for redshift, z"""
        
        H_0 = 100*self.h        #in km/s/Mpc    
    
        integral  = quad(lambda x: (self.c)/(H_0*numpy.sqrt(self.omega_m0*(1+x)**3 + 
            self.omega_lam0 + self.omega_r0*(1+z)**4 + (1 - self.omega0)*(1+x)**2)), 0, z)
        
        # flat universe
        if self.kappa == 0.0:
            return integral[0]
        # open universe
        elif (self.kappa < 0.0):
            return abs(self.kappa)**(-0.5) * numpy.sinh(abs(self.kappa)**(0.5)*integral[0])
        # closed universe
        elif (self.kappa > 0.0):
            return abs(self.kappa)**(-0.5) * numpy.sin(abs(self.kappa)**(0.5)*integral[0])
            
    def angDiameterDist(self, z):
        """computes angular diameter distance in Mpc for redshift, z"""        
        
        comoving = self.comovingDist(z)
                      
        angDists = comoving / (1+z)

        return angDists
        
    def luminosityDist(self, z):
        """computes luminosity distance in Mpc for redshift, z"""        
        
        comoving = self.comovingDist(z)
                      
        lumDists = comoving * (1+z)

        return lumDists


    def comovingVolume(self, z1, z2):
        """ computes comoving volume between redshifts z1 and z2 """

        H_0 = 100*self.h

        vol = quad(lambda x: (self.c/H_0)*(1+x)**2 * self.angDiameterDist(x)**2 / self.E(x),  z1, z2)
        
        return 4*numpy.pi*vol[0]
        
    def angularSize(self, z, diameter):
        """computes the angular size of an object of diameter d (in Mpc) and at redshift z"""
        
        angDists = self.angDiameterDist(z)
       
        angSize = diameter/angDists
        
        return angSize

    def magnitudeToLuminosity(self, Mr, kcorr):

        h = 0.71
        if (kcorr == 0.1):
            Mr_sun = 4.76 
        if (kcorr == 0.25):
            Mr_sun = 4.5
        if (kcorr == None):
            Mr_sun = 4.64

        lum = 10**((Mr_sun - Mr)/2.50)
        l = lum / 1e10 * h**2
    
        return l

    def getAgeAtZ(self, z):
        """
        return age of universe at redshift z in years
        """
    
        H_0 = 100*self.h/(978.46e9)     #in 1/years    
    
        a = 1./(1.+z)
        integral  = quad(lambda x: 1.0/(H_0*x*numpy.sqrt(self.omega_m0/x**3 + self.omega_lam0 
           + self.omega_r0/x**4 + (1 - self.omega0)/x**2)), 0, a)
            
        return integral[0] # in years
        
    def getAgeOfUniverse(self):
    
        return self.getAgeAtZ(0.0)