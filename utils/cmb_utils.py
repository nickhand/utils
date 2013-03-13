import numpy as np
from physical_constants import *

def planck_law(nu, T = T_cmb):
    """
    @brief Planck's law
    @param nu: frequency in GHz
    @keyword T: the temperature to use in Kelvin
    @return specific intensity in MJy/sr
    """
    # the dimensionaless frequency parameter
    x = h_planck*giga*nu / (k_b*T)
    
    # intensity in erg/s/cm/cm/Hz
    b_nu = (2*h_planck*((giga*nu)**3) / c_light**2) * 1/(np.exp(x)- 1.)
    
    # return in MJy/sr
    return  b_nu / (mega*jansky)


def rayleigh_jeans_law(nu, T = T_cmb):
    """
    @brief Rayleigh Jeans limit of Planck's law for h*nu << kT
    @param nu: frequency in GHz
    @keyword T: the temperature to use in Kelvin
    @return intensity in MJy/sr
    """
    
    # intensity in erg/s/cm/cm/Hz
    b_nu = 2. * (giga*nu)**2 * k_b * T / c_light**2
    
    # return in MJy/sr
    return b_nu / (mega*jansky)


def dBdT(nu, T = T_cmb):
    """
    @brief derivative of Planck's law with respect to temperature
    @param nu: frequency in GHz
    @keyword T: the temperature to use in Kelvin
    @return value of derivative in cgs units
    """
    # the constant out front
    A = 2 * h_planck**2 * (giga*nu)**4 / k_b  / (c_light*T)**2
    
    # dimensionaless freq param
    x = h_planck*giga*nu / k_b / T
    
    # return the derivative
    return  A * np.exp(x) / (np.exp(x) - 1.)**2
    

def dTdB(nu, T = T_cmb):
    """
    @brief inverse of the derivative of Planck's law with respect to temperature
    @param nu: frequency in GHz
    @keyword T: the temperature to use in Kelvin
    @return value of derivative in cgs units
    """

    return dB_dT(nu, T)**(-1.)

def dBdT_RJ(nu):
    """
    @brief derivative of the Rayleigh Jeans approximation with respect to temperature
    @param nu: frequency in GHz
    @keyword T: the temperature to use in Kelvin
    @return value of derivative in cgs units
    """

    return 2. * k_b * (giga*nu)**2 / c_light**2


def centroid(mp):
    """
    @brief compute the centroid of an input map, weighted by the data
    @param mp: the map object (flipper.flipper.liteMap)
    @return value of centroid in degrees
    """
    
    data = mp.data
    Nx = mp.Nx
    Ny = mp.Ny

    # compute the x pixels, weighted by data
    xsum = 0.0
    for i in xrange(0, Nx):
        for j in xrange(0, Ny):
            xsum += data[j, i] * (i+1)

    # compute the y pixels, weighted by data
    ysum = 0.0
    for i in range(0, Nx):
        for j in range(0, Ny):
            ysum += data[j, i] * (j+1)
      
    # the total weight normalization  
    totalsum = np.sum(data[:,:])
    
    # get the x, y weighted pixels
    xcm = xsum / totalsum
    ycm = ysum / totalsum
    
    # convert pixels to degrees
    cmdeg = mp.pixToSky(xcm, ycm)

    return cmdeg


def f_sz(nu):
    """
    @brief the frequency dependence of the thermal SZ effect
    @param nu: frequency in GHz
    @return value of frequency function
    """

    x = h_planck*giga*nu / k_b / T_cmb
    f = x*(np.exp(x) + 1.) / (np.exp(x) - 1.) - 4.0

    return f