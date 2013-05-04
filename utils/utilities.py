import numpy as np
import progressbar as pb
from scipy.special import gammaincc
from scipy.optimize import bisect
from scipy.integrate import quad
import pylab
from catIO import catalog
import collections
import sys
import itertools
import operator


def update_dict(d, value, keysToUpdate):
    """
    @brief update keys in a dictionary by strip replacing a value in the input dictionary
    """ 
    
    newDict = d.copy() # don't overwrite the input copy
    
    # update each key-value pair
    for k, v in keysToUpdate.iteritems():
        
        cnt = v.count("%s")
        replace = (value,)*cnt
        
        if k in newDict.keys():
            newDict[k] = v %replace    
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
                    x[key][int(fields[-1])] = v %replace
                else:
                    x[key] = v %replace     
                
            else:
                raise KeyError
        
    return newDict

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
    
    return
    
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
        
def initializeProgressBar(N, fd=sys.stderr):
    """
    @brief initialize a progress bar with N total ticks
    """
    bar = pb.ProgressBar(maxval=N, fd=fd, term_width = 100, widgets=[pb.Bar('=', '[', ']'), ' ', pb.Percentage()])
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
    
def smooth(x, window_len=10, window='hanning'):
    """
    @brief smooth the data using a window with requested size.

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

def gauss_kern(size, sizey=None):
    """ 
    @brief returns a normalized 2D gauss kernel array for convolutions 
    """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
    g = np.exp(-(x**2/float(size) + y**2/float(sizey)))
    return g / g.sum()

def blur_image(im, n, ny=None) :
    """ 
    @brief blurs the image by convolving with a gaussian kernel of typical
    size n. The optional keyword argument ny allows for a different
    size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve(im, g, mode='valid')
    return(improc)
    
def get_consecutive_ints(array):
    """
    @brief break the input array into arrays of consecutive integers
    """
    output = []
    for k, g in itertools.groupby(enumerate(array), lambda (i,x):i-x):
        output.append(map(operator.itemgetter(1), g))
        
    return np.array(output)
    
def chunks(l, n):
    """ 
    @brief yield n successive chunks from l.
    """
    out = []
    newn = int(len(l) / n)
    for i in xrange(0, n-1):
        out.append(l[i*newn:i*newn+newn])
    
    out.append(l[n*newn-newn:])
    
    return out
    
def latex_float(f):
    """
    @brief format a floating point number into a raw latex string
    """
    float_str = "{0:.2g}".format(f)
    if "e" in float_str:
        base, exponent = float_str.split("e")
        return r"{0} \ \times \ 10^{{{1}}}".format(base, int(exponent))
    else:
        return float_str
        


