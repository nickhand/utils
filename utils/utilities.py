import numpy as np
from physical_constants import *
import progressbar as pb
from scipy.special import gammaincc
from scipy.optimize import bisect
from scipy.integrate import quad
import pylab
from catIO import catalog
import collections
import sys

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


def bin(arrayX, arrayY, nBins, log = False, Nconst=False, norm=None, operator=np.mean):
    """
    @brief bin the inumpyut arrays using the given operator
    @param arrayX x input array
    @param arrayY y input array to be binned
    @param nBins number of bins to use
    @param log boolean denoting whether to use even bins in logspace
    @param Nconst boolean denoting whether to have fixed number of points per bin
    @param norm possible weights array to normalize binning by
    @param operator the operator to use when binning
    
    @return X binned X array
    @return Y binned Y array
    @return Ystd 1sigma errors on each Y bin
    @return weights number of points in each bin
    """
    
    # convert to arrays
    arrayX = np.array(arrayX)
    arrayY = np.array(arrayY)

    if Nconst == True:
        
        # define min and max distance and number of bins
        Nwidth = int(len(arrayY) / nBins)
        
    
        # initialize lists for later use 
        X = []          
        Y = []       
        Ystd = []
        weights = []     
        

        # sort arrays
        index = np.argsort(arrayX)
        arrayY = arrayY[index]
        arrayX = np.sort(arrayX)
        
    

        # create bins and calculate list values
        for i in range(0, nBins):
            
            x = arrayX[i*Nwidth:(i+1)*Nwidth]
            y = arrayY[i*Nwidth:(i+1)*Nwidth]

            if norm is not None:
                w = norm[i*Nwidth:(i+1)*Nwidth]
            
            X.append(np.mean(x))
            if norm is None:
                Y.append(operator(y))
            else:
                Y.append(np.sum(y)/np.sum(w))
                
            weights.append(len(x))
            
            Ystd.append(np.std(y))

    
        # converts lists to arrays
        X = np.array(X)
        Y = np.array(Y)
        Ystd = np.array(Ystd)
        weights = np.array(weights)

        
        return X, Y, Ystd, weights

    else:
        # define min and max distance and number of bins
        binWidth = (np.amax(arrayX) - np.amin(arrayX))/nBins

        max = np.amax(arrayX)
        min = np.amin(arrayX)

        if log:
            bins = np.logspace(np.log10(min), np.log10(max), nBins+1)
            
        # initialize lists for later use 
        X = []          
        Y = []       
        Ystd = []
        weights = []     
        

        # sort arrays
        index = np.argsort(arrayX)
        arrayY = arrayY[index]
        arrayX = np.sort(arrayX)

    
        # create bins and calculate list values
        for i in range(0, nBins):

            if log:
                if (i == nBins-1):
                    cond = np.where(arrayX >= bins[i])
                else:
                    cond = np.where((arrayX >= bins[i])*(arrayX < bins[i+1]))
            else:
                cut_low = min + i*binWidth
                cut_high = min + (i+1)*binWidth
        
                cond = np.where((arrayX >= cut_low)*(arrayX < cut_high))
            
            assert(len(cond[0]) > 0)
            
            x = arrayX[cond]
            y = arrayY[cond]
            if norm is not None:
                w = norm[cond]
            
            X.append(np.mean(x))
            if norm is None:
                Y.append(operator(y))
            else:
                Y.append(np.sum(y)/np.sum(w))
    
            weights.append(len(x))
            Ystd.append(np.std(y))

    
        # converts lists to arrays
        X = np.array(X)
        Y = np.array(Y)
        Ystd = np.array(Ystd)
        weights = np.array(weights)

        
        return X, Y, Ystd, weights

def extrap1d(interpolator):
    """
    @brief 1d extrapolator function, using linear extrapolation
    @param interpolator the interpolator function (ie, from scipy.interp1d)
    
    @return ufunclike the extrapolator function
    """
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0] + (x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1] + (x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return np.array(map(pointwise, np.array(xs)))

    return ufunclike

                                                                            
    
def getSigmaFromChiSquared(chi_sq, dof):
    """
    @brief compute the significance (in sigma) of a chi squared value
           and degrees of freedom (from Numerical Recipes)
    @param chi_sq the chi squared value
    @param dof the degrees of freedom
    
    @return the significance in sigma
    """
    def func1(x):
        return x - 1.0 + gammaincc(dof/2.0, chi_sq/2.0)
    
    # first calculate p, the confidence limit
    p = bisect(func1, 0, 1.0)

    def func2(x):

        val = quad(lambda y: 1.0/np.sqrt(2.0*np.pi)*np.exp(-y**2/2.0), -x, x)

        return val[0] - p

    # now calculate sigma
    sigma = bisect(func2, 0, 100)

    # 1 - p is probability that random variates could have this chi-squared value
    return sigma, 1-p
    
def getSigmaFromPValue(p_value):
    """
    @brief compute the significance (in sigma) of a given p-value
    @param p_value: 1 - p_value is prob that random variates could have a given chi squared

    @return the significance in sigma
    """

    def func(x):

        val = quad(lambda y: 1.0/np.sqrt(2.0*np.pi)*np.exp(-y**2/2.0), -x, x)

        return val[0] - p_value

    # now calculate sigma
    sigma = bisect(func, 0, 100)

    return sigma

def paper_single():
    """
    @brief initialize parameters for a single column pylab plot
    """
    
    pylab.rc('figure', figsize=(3.375,3.375))
    pylab.rc('figure.subplot', left=0.18, right=0.93, bottom=0.15, top=0.95)
    pylab.rc('lines', linewidth=1.5)
    pylab.rc('font', size=10.0)
    pylab.rc('xtick', labelsize='small')
    pylab.rc('ytick', labelsize='small')
    pylab.rc('legend', fontsize='medium') 
    pylab.rc('axes', linewidth=1.5)
    
def getIDFromRADec(ra, dec, tag):
    """
    @brief compute the ID in IAU format from ra, dec in decimal format
    """
    ra_s, dec_s = catalog.convertRADecDegreesToSexagesimal(ra, dec)
    tup = ra_s + dec_s

    if (dec < 0):
        iau_name = tag+"_J%02d%02d%05.2f-%02i%02d%05.1f" %tup
    else:
        iau_name = tag+"_J%02d%02d%05.2f+%02i%02d%05.1f" %tup

    return iau_name
        
def initializeProgressBar(N):
    """
    @brief initialize a progress bar with N total ticks
    """
    bar = pb.ProgressBar(maxval=N, term_width = 100, widgets=[pb.Bar('=', '[', ']'), ' ', pb.Percentage()])
    return bar
    
def getDuplicatesFromList(L):
    """
    @brief return the duplicate objects and indices of those objects from list L
    """
    dups = collections.defaultdict(list)
    for i, e in enumerate(L):
        dups[e].append(i)
        
    out = {}
    for k, v in sorted(dups.iteritems()):
        if len(v) >= 2:
            out[k] = v
    
    return out
    
def weighted_percentile(data, wt, percentiles): 
    """
    Compute weighted percentiles. 
    If the weights are equal, this is the same as normal percentiles. 
    Elements of the C{data} and C{wt} arrays correspond to 
    each other and must have equal length (unless C{wt} is C{None}). 
    
    @param data: The data. 
    @type data: A L{np.ndarray} array or a C{list} of numbers. 
    @param wt: How important is a given piece of data. 
    @type wt: C{None} or a L{np.ndarray} array or a C{list} of numbers. 
            All the weights must be non-negative and the sum must be 
            greater than zero. 
    @param percentiles: what percentiles to use.  (Not really percentiles, 
                        as the range is 0-1 rather than 0-100.) 
    @type percentiles: a C{list} of numbers between 0 and 1. 
    @rtype: [ C{float}, ... ] 
    @return: the weighted percentiles of the data. 
    """ 
    assert np.greater_equal(percentiles, 0.0).all(), "Percentiles less than zero" 
    assert np.less_equal(percentiles, 1.0).all(), "Percentiles greater than one" 
    data = np.asarray(data) 
    assert len(data.shape) == 1 
    if wt is None: 
        wt = np.ones(data.shape, np.float) 
    else: 
        wt = np.asarray(wt, np.float) 
        assert wt.shape == data.shape 
        assert np.greater_equal(wt, 0.0).all(), "Not all weights are non-negative." 
    assert len(wt.shape) == 1 
    n = data.shape[0] 
    assert n > 0 
    i = np.argsort(data) 
    sd = np.take(data, i, axis=0) 
    sw = np.take(wt, i, axis=0) 
    aw = np.add.accumulate(sw) 
    if not aw[-1] > 0: 
        raise ValueError, "Nonpositive weight sum" 
    w = (aw-0.5*sw)/aw[-1] 
    spots = np.searchsorted(w, percentiles) 
    o = [] 
    for (s, p) in zip(spots, percentiles): 
        if s == 0: 
            o.append(sd[0]) 
        elif s == n: 
            o.append(sd[n-1]) 
        else: 
            f1 = (w[s] - p)/(w[s] - w[s-1]) 
            f2 = (p - w[s-1])/(w[s] - w[s-1]) 
            assert f1>=0 and f2>=0 and f1<=1 and f2<=1 
            assert abs(f1+f2-1.0) < 1e-6 
            o.append(sd[s-1]*f1 + sd[s]*f2) 
    return o 
    
def stringToFunction(astr):
    """
    @brief given a string containing the name of a function, convert 
    it to a function

    @param astr: string of function name
    """
    module, _, function = astr.rpartition('.')
    if module:
        __import__(module)
        mod = sys.modules[module]
    else:
        mod = sys.modules['__main__']

    return getattr(mod, function)
    
def runge_kutta_4th(x, y, z, a , dx):
    """
    @brief a fourth order Runge-Kutta ODE solver
    @param x: the independent variable (float)
    @param y: the initial dependent variable at x (float)
    @param z: the first derivative dy/dx (float)
    @param a: the acceleration function of y, d2y/dx2, a(x, y, z) (function)
    @param dx: the interval step in x (float)
    """

    x1 = x
    y1 = y
    z1 = z
    a1 = a(x1, y1, z1) #this is dz/dx

    x2 = x + 0.5*dx
    y2 = y + 0.5*z1*dx
    z2 = z + 0.5*a1*dx
    a2 = a(x2, y2, z2)

    x3 = x + 0.5*dx
    y3 = y + 0.5*z2*dx
    z3 = z + 0.5*a2*dx
    a3 = a(x3, y3, z3)

    x4 = x + dx
    y4 = y + z3*dx
    z4 = z + a3*dx
    a4 = a(x3, y4, z4)

    yf = y + (dx/6.0)*(z1 + 2.0*z2 + 2.0*z3 + z4)
    zf = z + (dx/6.0)*(a1 + 2.0*a2 + 2.0*a3 + a4)
    xf = z + dx

    return yf, zf
    
def weighted_mean_arrays(vals, errs):
    """
    @brief read in a list of lists and errors and compute the mean
    """
    vals = np.asarray(vals)
    errs = np.asarray(errs)
    
    mean_vals = np.sum(vals/errs**2, axis=0) / np.sum(1/errs**2, axis=0)
    mean_errs = np.sum(1./errs**2, axis=0)**(-0.5)
    
    return mean_vals, mean_errs
    
