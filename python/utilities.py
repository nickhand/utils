import numpy
import scipy
import random
from scipy.special import gammaincc
from scipy.optimize import bisect
from scipy.integrate import quad
import pylab
from physical_constants import *
from flipper import *
import catalog
import progressbar as pb
import collections


def B_nu(nu, T = T_CMB):
    """
    @brief Planck's law
    @param nu frequency in GHz
    @return intensity in MJy/sr
    """
   
    x = h_planck*nu*1e9 / (k_b*T)
    b_nu = (2*h_planck*((nu*1e9)**3) / c_light**2) * 1/(numpy.exp(x)- 1.)
    
    return  b_nu * cgsToJansky / 1e6



def RJ(nu, T = T_CMB):
    """
    @brief Rayleigh Jeans limit of Planck's law for h*nu << kT
    @param nu frequency in GHz
    @return intensity in MJy/sr
    """
    b_nu = 2. * (nu*1e9)**2 * k_b * T / c_light**2
    return b_nu * cgsToJansky / 1e6


def dB_dT(nu, T = T_CMB):
    """
    @brief Derivative of Planck's Law with respect to temperature
    @param nu frequency in GHz
    @return value of derivative in cgs units
    """
    A = 2 * h_planck**2 * (nu*1e9)**4 / k_b  / (c_light*T)**2
    
    x = h_planck*nu*1e9 / k_b / T
    dB_dT = A * numpy.exp(x) / (numpy.exp(x) - 1.)**2
    
    return dB_dT

def dT_dB(nu, T = T_CMB):
    """
    @brief inverse of the derivative of Planck's Law with respect to temperature
    @param nu frequency in GHz
    @return value of derivative in cgs units
    """

    return dB_dT(nu, T)**(-1.0)


def dB_dTrj(nu):
    """
    @brief Derivative of the RJ approximation with respect to temperature
    @param nu frequency in GHz
    @return value of derivative in cgs units
    """

    dB_dTrj = 2. * k_b * (nu*1e9)**2 / c_light**2

    return dB_dTrj



def centroid(mp):
    """
    @brief calculates the centroid of map
    @param mp liteMap object
    @return value of centroid in degrees
    """
    
    data = mp.data
    Nx = mp.Nx
    Ny = mp.Ny

    xsum = 0.0
    for i in numpy.range(0, Ny):
        for j in numpy.range(0, Nx):
            xsum += data[i, j] * (j+1)

    ysum = 0.0
    for i in range(0, Nx):
        for j in range(0, Ny):
            ysum += data[j, i] * (j+1)
        
    totalsum = numpy.sum(data[:,:])
    xcm = xsum / totalsum
    ycm = ysum / totalsum
    cmdeg = mp.pixToSky(xcm, ycm)

    return cmdeg


def f_sz (nu):
    """
    @brief calculates frequency dependence of thermal SZ effect
    @param nu freq in GHz
    @return value of frequency function
    """

    x = h_planck*nu*1e9 / k_b / T_CMB
    f = x*(numpy.exp(x) + 1.) / (numpy.exp(x) - 1.) - 4.0

    return f

def h_sz (nu):

    x = h_planck*nu*1e9 / k_b / T_cmb
    g =  f_sz(nu) * x * numpy.exp(x) / (numpy.exp(x) - 1.)

    return g

    
def convolveWithBeam(mp, beamFile, template = None, pixelPitch = None):
    """
    @brief convolve litemap object with ACT beam
    @param mp inumpyut litemap object
    @param beamFile filename of the beam data
    @param template template litemap object
    @param pixelPitch pixelPitch of mp
    @return convolved map object
    """

    
    f = file(beamFile)
    ell=[]
    wl = []
    for line in f:
        fields = line[:-1].split()
        ell.append(float(fields[0]))
        wl.append(float(fields[1]))

    ell = numpy.array(ell)
    wl  = numpy.array(wl)

    if (template == None):
        t = mp.copy()
    else:
        t = template.copy()
    t.data[:] = 0.0
    ft = fftTools.fftFromLiteMap(t)

    l_f = numpy.floor(ft.modLMap)
    l_c = numpy.ceil(ft.modLMap)

    for i in xrange(numpy.shape(ft.kMap)[0]):
        for j in xrange(numpy.shape(ft.kMap)[1]):
            w_lo = wl[l_f[i,j]]
            w_hi = wl[l_c[i,j]]
            trueL = ft.modLMap[i,j]
            w = (w_hi-w_lo)*(trueL - l_f[i,j]) + w_lo
            ft.kMap[i,j] = w

    t.data = numpy.sqrt(abs(ft.kMap))

    mp_ft = fftTools.fftFromLiteMap(mp)

    if (pixelPitch != None):
        t = liteMap.upgradePixelPitch(t, pixelPitch)

    mp_ft.kMap *= t.data
    data_conv = numpy.real(numpy.fft.ifft2(mp_ft.kMap))
    t.data[:] = data_conv[:]

    return t

def mask(lm, ra, dec, hw = 7):
    """
    @brief read in a litemap and implement a mask at the given ra/dec
    @param lm the litemap to mask
    @param ra the right ascension of the mask location
    @param dec the declination of the mask location
    @param hw the halfwidth of the mask in pixels
    
    @return masked the masked litemap
    """


    x = range(0, lm.Nx)
    y = range(0, lm.Ny)

    xx, yy = numpy.meshgrid(x, y)

    x0, y0 = lm.skyToPix(ra, dec)
    dist = numpy.sqrt((xx - x0)**2 + (yy - y0)**2)

    a,b = numpy.where((dist > 15) * (dist < 25))
    avg = numpy.mean(lm.data[a,b])

    i, j = numpy.where(dist < hw)
    lm.data[i,j] = avg

    masked = liteMap.liteMapFromDataAndWCS(lm.data, lm.wcs)

    return masked


def bin(arrayX, arrayY, nBins, log = False, Nconst=False, norm=None, operator=numpy.mean):
    """
    @brief bin the inumpyut arrays using the given operator
    @param arrayX x inumpyut array
    @param arrayY y inumpyut array to be binned
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
    arrayX = numpy.array(arrayX)
    arrayY = numpy.array(arrayY)

    if Nconst == True:
        
        # define min and max distance and number of bins
        Nwidth = int(len(arrayY) / nBins)
        
    
        # initialize lists for later use 
        X = []          
        Y = []       
        Ystd = []
        weights = []     
        

        # sort arrays
        index = numpy.argsort(arrayX)
        arrayY = arrayY[index]
        arrayX = numpy.sort(arrayX)
        
    

        # create bins and calculate list values
        for i in range(0, nBins):
            
            x = arrayX[i*Nwidth:(i+1)*Nwidth]
            y = arrayY[i*Nwidth:(i+1)*Nwidth]

            if norm is not None:
                w = norm[i*Nwidth:(i+1)*Nwidth]
            
            X.append(numpy.mean(x))
            if norm is None:
                Y.append(operator(y))
            else:
                Y.append(numpy.sum(y)/numpy.sum(w))
                
            weights.append(len(x))
            
            Ystd.append(numpy.std(y))

    
        # converts lists to arrays
        X = numpy.array(X)
        Y = numpy.array(Y)
        Ystd = numpy.array(Ystd)
        weights = numpy.array(weights)

        
        return X, Y, Ystd, weights

    else:
        # define min and max distance and number of bins
        binWidth = (numpy.amax(arrayX) - numpy.amin(arrayX))/nBins

        max = numpy.amax(arrayX)
        min = numpy.amin(arrayX)

        if log:
            bins = numpy.logspace(numpy.log10(min), numpy.log10(max), nBins+1)
            
        # initialize lists for later use 
        X = []          
        Y = []       
        Ystd = []
        weights = []     
        

        # sort arrays
        index = numpy.argsort(arrayX)
        arrayY = arrayY[index]
        arrayX = numpy.sort(arrayX)

    
        # create bins and calculate list values
        for i in range(0, nBins):

            if log:
                if (i == nBins-1):
                    cond = numpy.where(arrayX >= bins[i])
                else:
                    cond = numpy.where((arrayX >= bins[i])*(arrayX < bins[i+1]))
            else:
                cut_low = min + i*binWidth
                cut_high = min + (i+1)*binWidth
        
                cond = numpy.where((arrayX >= cut_low)*(arrayX < cut_high))
            
            assert(len(cond[0]) > 0)
            
            x = arrayX[cond]
            y = arrayY[cond]
    
        
            X.append(numpy.mean(x))
            Y.append(operator(y))
            weights.append(len(x))
            Ystd.append(numpy.std(y))

    
        # converts lists to arrays
        X = numpy.array(X)
        Y = numpy.array(Y)
        Ystd = numpy.array(Ystd)
        weights = numpy.array(weights)

        
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
        return numpy.array(map(pointwise, numpy.array(xs)))

    return ufunclike


def getBootstrapError(x, w,  N_resamples):
    """
    @brief calculates the errors through bootstrapping
    @param x array to be resampled from
    @param w possible weight array
    @param N_resamples number of resamples
    
    @return mean of resampled distribution
    @return 1sigma error on mean of resampled dist.
    """

    mean_list = []

    N_bin = len(x)

    if w is None:
        w = numpy.ones(N_bin)
        
    for j in xrange(0, N_resamples):

        avg = 0
        norm = 0
        for k in xrange(0, N_bin):
            
            i = random.randint(0, N_bin-1)
            avg += x[i]*w[i]
            norm += abs(w[i])
           
        
        mean_list.append(avg/norm)
            
    return numpy.mean(mean_list), numpy.std(mean_list)
                                                                            
    
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

        val = quad(lambda y: 1.0/numpy.sqrt(2.0*numpy.pi)*numpy.exp(-y**2/2.0), -x, x)

        return val[0] - p

    # now calculate sigma
    sigma = bisect(func2, 0, 100)

    # 1 - p is probability that random variates could have this chi-squared value
    return sigma, 1-p

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
    @brief print out the duplicate objects and indices of those objects from list L
    """
    dups = collections.defaultdict(list)
    for i, e in enumerate(L):
        dups[e].append(i)
        
    for k, v in sorted(dups.iteritems()):
        if len(v) >= 2:
            print '%s: %r' %(k, v)
    
    return 
    
def weighted_percentile(data, wt, percentiles): 
    """
    Compute weighted percentiles. 
    If the weights are equal, this is the same as normal percentiles. 
    Elements of the C{data} and C{wt} arrays correspond to 
    each other and must have equal length (unless C{wt} is C{None}). 
    
    @param data: The data. 
    @type data: A L{numpy.ndarray} array or a C{list} of numbers. 
    @param wt: How important is a given piece of data. 
    @type wt: C{None} or a L{numpy.ndarray} array or a C{list} of numbers. 
            All the weights must be non-negative and the sum must be 
            greater than zero. 
    @param percentiles: what percentiles to use.  (Not really percentiles, 
                        as the range is 0-1 rather than 0-100.) 
    @type percentiles: a C{list} of numbers between 0 and 1. 
    @rtype: [ C{float}, ... ] 
    @return: the weighted percentiles of the data. 
    """ 
    assert numpy.greater_equal(percentiles, 0.0).all(), "Percentiles less than zero" 
    assert numpy.less_equal(percentiles, 1.0).all(), "Percentiles greater than one" 
    data = numpy.asarray(data) 
    assert len(data.shape) == 1 
    if wt is None: 
        wt = numpy.ones(data.shape, numpy.float) 
    else: 
        wt = numpy.asarray(wt, numpy.float) 
        assert wt.shape == data.shape 
        assert numpy.greater_equal(wt, 0.0).all(), "Not all weights are non-negative." 
    assert len(wt.shape) == 1 
    n = data.shape[0] 
    assert n > 0 
    i = numpy.argsort(data) 
    sd = numpy.take(data, i, axis=0) 
    sw = numpy.take(wt, i, axis=0) 
    aw = numpy.add.accumulate(sw) 
    if not aw[-1] > 0: 
        raise ValueError, "Nonpositive weight sum" 
    w = (aw-0.5*sw)/aw[-1] 
    spots = numpy.searchsorted(w, percentiles) 
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
    
