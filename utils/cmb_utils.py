"""
 cmb_utils.py
 utility functions useful for CMB analysis
 
 author: Nick Hand
 contact: nhand@berkeley.edu
"""
import numpy as np
from physical_constants import *

def planck_law(nu, T=T_cmb):
    """
    Planck's law, B_nu
    
    Parameters
    ----------
    nu : float or numpy.ndarray 
        the frequency in GHz
    T : float, optional 
        the temperature to evaluate Planck's law at in Kelvin. Default is T_cmb.
        
    Returns
    -------
    B_nu : float or numpy.ndarray
        specific intensity in MJy/sr
    """
    x = h_planck*giga*nu / (k_b*T)
    b_nu = (2*h_planck*((giga*nu)**3) / c_light**2) * 1/(np.exp(x)- 1.)
    return  b_nu / (mega*jansky)
#end planck_law

#-------------------------------------------------------------------------------
def rayleigh_jeans_law(nu, T=T_cmb):
    """
    The Rayleigh Jeans limit of Planck's law for h*nu << kT
    
    Parameters
    ----------
    nu : float or numpy.ndarray 
        the frequency in GHz
    T : float, optional 
        the temperature to evaluate Planck's law at in Kelvin. Default is T_cmb.
        
    Returns
    -------
    B_nu : float or numpy.ndarray
        specific intensity in MJy/sr
    """
    b_nu = 2. * (giga*nu)**2 * k_b * T / c_light**2
    return b_nu / (mega*jansky)
#end rayleigh_jeans_law

#-------------------------------------------------------------------------------
def dBdT(nu, T=T_cmb):
    """
    The derivative of Planck's law with respect to temperature.
    
    Parameters
    ----------
    nu : float or numpy.ndarray 
        the frequency in GHz
    T : float, optional 
        the temperature to evaluate Planck's law at in Kelvin. Default is T_cmb.
    
    Returns
    -------
    deriv : float or numpy.ndarray
        the value of derivative in cgs units
    """
    A = 2 * h_planck**2 * (giga*nu)**4 / k_b  / (c_light*T)**2
    x = h_planck*giga*nu / k_b / T
    return  A * np.exp(x) / (np.exp(x) - 1.)**2
#end dBdT

#-------------------------------------------------------------------------------
def dTdB(nu, T = T_cmb):
    """
    The inverse of the derivative of Planck's law with respect to temperature.
    
    Parameters
    ----------
    nu : float or numpy.ndarray 
        the frequency in GHz
    T : float, optional 
        the temperature to evaluate Planck's law at in Kelvin. Default is T_cmb.
    
    Returns
    -------
    deriv : float or numpy.ndarray
        the value of derivative in cgs units
    """

    return dB_dT(nu, T)**(-1.)
#end dTdB

#-------------------------------------------------------------------------------
def dBdT_RJ(nu):
    """
    The derivative of the Rayleigh Jeans approximation with respect to 
    temperature, which is independent of temperature.
    
    Parameters
    ----------
    nu : float or numpy.ndarray 
        the frequency in GHz
    T : float, optional 
        the temperature to evaluate Planck's law at in Kelvin. Default is T_cmb.
    
    Returns
    -------
    deriv : float or numpy.ndarray
        the value of derivative in cgs units
    """
    return 2. * k_b * (giga*nu)**2 / c_light**2
#end dBdT_RJ

#-------------------------------------------------------------------------------
def centroid(mp):
    """
    Compute the centroid of an input map, weighted by the data
    
    Parameters
    ----------
    mp : flipper.liteMap
        the input liteMap object storing the data
    
    Returns
    -------
    cmdeg : tuple
        the value of the centroid (ra, dec) in degrees
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
#end centroid

#-------------------------------------------------------------------------------
def f_sz(nu):
    """
    The frequency dependence of the thermal SZ effect
    
    Parameters
    ----------
    nu : float or numpy.ndarray 
        the frequency in GHz
    """
    x = h_planck*giga*nu / k_b / T_cmb
    return x*(np.exp(x) + 1.) / (np.exp(x) - 1.) - 4.0
#end f_sz

#-------------------------------------------------------------------------------