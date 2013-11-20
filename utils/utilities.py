import numpy as np
import progressbar as pb
import pylab
from catIO import catalog
import sys
import itertools
import operator
from utils import pytools
import collections

#-------------------------------------------------------------------------------
def vincenty_sphere_dist(ra1, dec1, ra2, dec2):
    """
    Method implementing the Vincenty formula for angular distance 
    on a sphere: stable at poles and antipodes but more 
    complex/computationally expensive. Input ra/dec are in degrees.

    
    Returns the angular separation in degrees.
    """
    toRad = np.pi/180.
    lon1, lat1 = ra1*toRad, dec1*toRad
    lon2, lat2 = ra2*toRad, dec2*toRad

    sdlon = np.sin(lon2 - lon1)
    cdlon = np.cos(lon2 - lon1)
    slat1 = np.sin(lat1)
    slat2 = np.sin(lat2)
    clat1 = np.cos(lat1)
    clat2 = np.cos(lat2)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denominator = slat1 * slat2 + clat1 * clat2 * cdlon

    return 180./np.pi*np.arctan2((num1**2 + num2**2)**0.5, denominator)
#end vincenty_sphere_dist

#-------------------------------------------------------------------------------
def update_dict(d, value, keysToUpdate):
    """
    Update keys in a dictionary by string replacing a value in the input dictionary
    
    Parameters
    ----------
    d : dict
        the dictionary to update
    value : int, str
        the value to insert into the dict
    keysToUpdate : dict
        dictionary with (k, v) where k is a key in d and v is the format
        string that value will be inserted into
    """ 
    newDict = d.copy()
    
    # update each key-value pair
    for k, v in keysToUpdate.iteritems():
        
        cnt = v.count("%s")
        replace = (value,)*cnt
        
        if k in newDict.keys():
            thistype = type(newDict[k])
            newDict[k] = thistype(v %replace)
        else:
            if '/' in k:
                fields = k.split('/')
                updated = False
                x = newDict
                index = 0
                if fields[-1].isdigit():
                    key = fields[-2]
                else:
                    key = fields[-1]
                while key not in x.keys():
                    x = x[fields[index]]
                    index += 1
            
                if fields[-1].isdigit():
                    thistype = type(x[key][int(fields[-1])])
                    x[key][int(fields[-1])] = thistype(v %replace)
                else:
                    thistype = type(x[key])
                    x[key] = thistype(v %replace)
                
            else:
                raise KeyError
        
    return newDict
#end update_dict

#-------------------------------------------------------------------------------
def bin(arrayX, arrayY, nBins, log = False, Nconst=False, norm=None, operator=np.mean):
    """
    Bin the input arrays using the given operator
    
    Parameters
    ----------
    arrayX : numpy.ndarray
        the x input array
    arrayY : numpy.ndarray 
        the input array to be binned
    nBins : int
        number of bins to use
    log : bool, optional 
        whether to use even bins in logspace or realspace
    Nconst : bool, optional
        whether to have fixed number of points per bin
    norm : numpy.ndarray
        a possible weights array to normalize binning by
    operator : function 
        the operator to use when binning
    
    Returns
    -------
    X : numpy.ndarray
        the binned X array
    Y : numpy.ndarray 
        the binned Y array
    Ystd : numpy.ndarray
        the 1-sigma errors on the value in each bin
    weights : numpy.ndarray
        the number of points in each bin
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
#end bin

#-------------------------------------------------------------------------------
def extrap1d(interpolator):
    """
    A 1d extrapolator function, using linear extrapolation
    
    Parameters
    ----------
    interpolator : scipy.interpolate.interp1d 
        the interpolator function
    
    Returns
    -------
    ufunclike : function
        the extrapolator function
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
#end extrap1d

#-------------------------------------------------------------------------------
def paper_single():
    """
    Initialize parameters for a single column pylab plot
    """
    
    pylab.rc('figure', figsize=(3.375,3.375))
    pylab.rc('figure.subplot', left=0.18, right=0.93, bottom=0.15, top=0.95)
    pylab.rc('lines', linewidth=1.5)
    pylab.rc('font', size=10.0)
    pylab.rc('xtick', labelsize='small')
    pylab.rc('ytick', labelsize='small')
    pylab.rc('legend', fontsize='medium') 
    pylab.rc('axes', linewidth=1.5)
#end paper_single
    
#-------------------------------------------------------------------------------
def getIDFromRADec(ra, dec, tag):
    """
    Compute the ID in IAU format from ra, dec in decimal format
    """
    ra_s, dec_s = catalog.convertRADecDegreesToSexagesimal(ra, dec)
    tup = ra_s + dec_s

    if (dec < 0):
        iau_name = tag+"_J%02d%02d%05.2f-%02i%02d%05.1f" %tup
    else:
        iau_name = tag+"_J%02d%02d%05.2f+%02i%02d%05.1f" %tup

    return iau_name
#end getIDFromRADec

#-------------------------------------------------------------------------------        
def initializeProgressBar(N, fd=sys.stderr):
    """
    Initialize a progress bar with N total ticks
    """
    return pb.ProgressBar(maxval=N, fd=fd, term_width=100, 
                         widgets=[pb.Bar('=', '[', ']'), ' ', pb.Percentage()])
#end initializeProgressBar
    
#-------------------------------------------------------------------------------
def getDuplicatesFromList(L):
    """
    Return the duplicate objects and indices of those objects from list L

    Parameters
    ----------
    L : list
        the list to find duplicates in
        
    Returns
    -------
    out : dict
        the dictionary with keys that are the duplicate items and values that
        are the list of indices of those items in L
    """
    dups = collections.defaultdict(list)
    for i, e in enumerate(L):
        dups[e].append(i)
        
    out = {}
    for k, v in sorted(dups.iteritems()):
        if len(v) >= 2:
            out[k] = v
    
    return out
#end getDuplicatesFromList
    
#-------------------------------------------------------------------------------
def stringToFunction(astr):
    """
    Given a string containing the name of a function, convert it to a function
    
    Parameters
    ----------
    astr : str 
        the string to convert to a function name
    """
    module, _, function = astr.rpartition('.')
    if module:
        __import__(module)
        mod = sys.modules[module]
    else:
        mod = sys.modules['__main__']

    return getattr(mod, function)
#end stringToFunction
    
#-------------------------------------------------------------------------------
def runge_kutta_4th(x, y, z, a , dx):
    """
    A fourth order Runge-Kutta ODE solver
    
    Parameters
    ----------
    x : float 
        the independent variable
    y : float
        the initial dependent variable at x
    z : the first derivative dy/dx
        float
    a : function
        the acceleration function of y, d2y/dx2, a(x, y, z)
    dx : float 
        the interval step in x
    
    Returns
    -------
    yf : float
        the updated value of y at x + dx
    zf : float
        the updated value of z at x + dx
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
#end runge_kutta_4th
    
#-------------------------------------------------------------------------------
def weighted_mean_arrays(vals, errs):
    """
    Read in lists of lists giving values and errors and compute weighted mean
    """
    vals = np.asarray(vals)
    errs = np.asarray(errs)
    
    inds = np.where(np.isnan(errs))
    errs[inds] = np.inf
    
    mean_vals = np.sum(vals/errs**2, axis=0) / np.sum(1/errs**2, axis=0)
    mean_errs = np.sum(1./errs**2, axis=0)**(-0.5)
    
    return mean_vals, mean_errs
#end weighted_mean_arrays

#-------------------------------------------------------------------------------    
def smooth(x, window_len=10, window='hanning'):
    """
    Smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal 
        window_len: the dimension of the smoothing window
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    import numpy as np    
    t = np.linspace(-2,2,0.1)
    x = np.sin(t)+np.random.randn(len(t))*0.1
    y = smooth(x)

    see also: 

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is one of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"

    s=np.r_[2*x[0]-x[window_len:1:-1], x, 2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))

    if window == 'flat': #moving average
        w = np.ones(window_len,'d')
    else:
        w = getattr(np, window)(window_len)
    y = np.convolve(w/w.sum(), s, mode='same')
    return y[window_len-1:-window_len+1]
#end smooth

#-------------------------------------------------------------------------------
def gauss_kern(size, sizey=None):
    """ 
    Returns a normalized 2D gauss kernel array for convolutions 
    """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
    g = np.exp(-(x**2/float(size) + y**2/float(sizey)))
    return g / g.sum()

#end gauss_kern

#-------------------------------------------------------------------------------
def blur_image(im, n, ny=None) :
    """ 
    Blurs the image by convolving with a gaussian kernel of typical
    size n. The optional keyword argument ny allows for a different
    size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im, g, mode='valid')
    return(improc)
#end blur_image

#-------------------------------------------------------------------------------
def get_consecutive_ints(array):
    """
    Break the input array into arrays of consecutive integers
    """
    output = []
    for k, g in itertools.groupby(enumerate(array), lambda (i,x):i-x):
        output.append(map(operator.itemgetter(1), g))
        
    return np.array(output)
#end get_consecutive_ints
    
#-------------------------------------------------------------------------------
def chunks(l, n):
    """ 
    Return n successive chunks from l.
    """
    out = []
    newn = int(len(l) / n)
    for i in xrange(0, n-1):
        out.append(l[i*newn:i*newn+newn])
    
    out.append(l[n*newn-newn:])
    
    return out
#end chunks
    
#-------------------------------------------------------------------------------
def latex_float(f):
    """
    Format a floating point number into a raw latex string
    """
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \ \times \ 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str
#end latex_float
        
#-------------------------------------------------------------------------------
def flatten(l):
    """
    A generator to flatten a list of irregular lists
    """
    for el in l:
        if isinstance(el, collections.Iterable) and not isinstance(el, basestring):
            for sub in flatten(el):
                yield sub
        else:
            yield el
#end flatten

#-------------------------------------------------------------------------------
def format_float(x, str_format, min_val=-3, max_val=3):
    """
    Format a float number
    """
    logx = np.log10(abs(x))
    if logx < min_val or logx > max_val:
        return (str_format+"e") %x
    else:
        return (str_format+"f") %x


