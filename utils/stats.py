"""
 stats.py
 utils: module with useful statistics functions
 
 author: Nick Hand
 contact: nhand@berkeley.edu
 creation date: 11/20/2013
"""
import numpy as np
import mpfit
import scipy.integrate as sint
import scipy.optimize as sopt
from scipy.interpolate import InterpolatedUnivariateSpline
import collections
from scipy.special import gammaincc, erf

def gaussian_likelihood(A, d, th, cov, log=False):
    """
    A gaussian likelihood function
    
    Parameters
    ----------
    A : float 
        the amplitude that multiplies the theory curve
    d : numpy.ndarray
        the data array
    th : numpy.ndarray
        the theory model
    cov : numpy.ndarray
        the covariance matrix
    log : bool, optional
        whether to return the log likelihood; default = False
    """
    C_inv = np.linalg.inv(cov) 
    Y_mu = d - A*th
    chisq = np.dot(Y_mu, np.dot(C_inv, Y_mu))
    det = np.linalg.det(cov)
    logp = -0.5*len(d)*np.log(2*np.pi*det)  - 0.5*chisq
    if log:
        return logp
    else:
        return np.exp(logp)
#end gaussian_likelihood
    
#-------------------------------------------------------------------------------
def amplitude_likelihood(A, data, model, covar):
    """
    Return the normalized likelihood given the input data, model, and covariance
    matrix for a series of model amplitudes. The likelihood is normalized
    such that its integral is unity
    
    Parameters
    ----------
    A : float or numpy.ndarray
        the parameter values to evaluate the likelihood at
    data : numpy.ndarray
        the data array
    model : numpy.ndarray
        the model array to fit the amplitude of
    covar : numpy.ndarray
        the covariance matrix specifying the uncertainties
    """
    if np.isscalar(A):
        return gaussian_likelihood(A, data, model, covar)/N
    else:
        p = np.array([gaussian_likelihood(x, data, model, covar) for x in A])
        N = sint.simps(p/p.max(), x=A)
        return p/p.max()/N
#end amplitude_likelihood
    
#-------------------------------------------------------------------------------
def fit_amplitude_ML(data, model, covar, quiet=False, compute_error=True):
    """
    Fit the data to a model using maximum-likelihood method, returning the
    best fit amplitude and error
    
    Parameters
    ----------
    data : numpy.ndarray
        the data array
    model : numpy.ndarray
        the model array to fit the amplitude of
    covar : numpy.ndarray
        the covariance matrix specifying the uncertainties
    quiet : bool, optional
        whether to quiet the convergence messages; default = False
    compute_error : bool, optional
        whether to compute the standard deviation; default = True
    
    Returns
    -------
    A : float
        the amplitude that maximized the likelihood
    A_err : float
        the standard deviation of the likelihood
    """        
    # do the minimization    
    quiet = int(not quiet)
    A_best = sopt.fmin(lambda x: -gaussian_likelihood(x, data, model, covar, True), 1.0, disp=quiet)[0]
    
    A_err = None
    if compute_error:
        
        # compute the proper normalization, such that the integral is unity
        p_best = gaussian_likelihood(A_best, data, model, covar)
        norm = sint.quad(lambda x: gaussian_likelihood(x, data, model, covar)/p_best, A_best, np.inf)[0]
        norm += sint.quad(lambda x: gaussian_likelihood(x, data, model, covar)/p_best, -np.inf, A_best)[0]
        norm *= p_best
        
        # compute the mean value
        mu = sint.quad(lambda x: x*gaussian_likelihood(x, data, model, covar)/norm, A_best, np.inf)[0]
        mu += sint.quad(lambda x: x*gaussian_likelihood(x, data, model, covar)/norm, -np.inf, A_best)[0]
      
        # compute the variances
        variance = sint.quad(lambda x: (x-mu)**2*gaussian_likelihood(x, data, model, covar)/norm, A_best, np.inf)[0]
        variance += sint.quad(lambda x: (x-mu)**2*gaussian_likelihood(x, data, model, covar)/norm, -np.inf, A_best)[0]
        A_err = np.sqrt(variance)
    
    return A_best, A_err
#end fit_amplitude_ML

#-------------------------------------------------------------------------------
def fit_amplitude(data, errs, model):
    """
    Fit an amplitude to the data using a given model, using a chi-squared
    least squares fit
    
    Parameters
    ----------
    data : numpy.ndarray
        the data array
    errs : numpy.ndarray
        the uncertainities on the data values
    model : numpy.ndarray
        the model array to fit the amplitude of
        
    Returns
    -------
    A : float
        the best-fit amplitude
    A_err : float
        the 1 sigma uncertainty in the amplitude
    """
    def modelFunction(theory, p):
        return p[0]*theory

    def fittingFunction(p, fjac=None, theory=None, y=None, err=None):
        model = modelFunction(theory, p)
        status = 0
        return ([status, (y-model)/err])

    # measure the amplitude
    p0 = [1.0]
    fa = {'theory': model, 'y': data, 'err': errs}
    m = mpfit.mpfit(fittingFunction, p0, functkw=fa)
    
    # the amplitude
    A = m.params[0]
    
    # degrees of freedom
    dof = len(data) - len(m.params) 
    
    # scaled uncertainties
    A_err = m.perror[0] * np.sqrt(m.fnorm / dof)
    
    return A, A_err
#end fit_amplitude

#-------------------------------------------------------------------------------
def random_variates_from_dist(x, dist, size=1, cumulative=False):
    """
    Returns a random variate from the input discrete distribution using the
    inverse transform sampling method
    
    Parameters
    ----------
    x : numpy.ndarray
        the random variates corresponding to dist
    dist : numpy.ndarray
        the desired probability distribution
    size : int, optional
        the number of random_variates to return; default = 1
    cumulative : bool, optional
        whether the input distribution is a cumulative distribution or not;
        default = False
    """
    if not cumulative:
        F = dist.cumsum() / sum(dist)
    else:
        F = 1.*dist/dist.max()
        
    f_interp = InterpolatedUnivariateSpline(x, F)
    xmin = np.amin(x)    
    xmax = np.amax(x)     

    if size == 1:    
        r = np.random.random(size=size)*(1.-xmin) + xmin
        return sopt.bisect(lambda a: f_interp(a)-r, xmin, xmax)
    else:
        rands = np.random.random(size=size)*(1.-xmin) + xmin
        ans = []
        for r in rands:
            ans.append(sopt.bisect(lambda a: f_interp(a)-r, xmin, xmax) )
        return np.array(ans)
#end random_variate_from_dist

#-------------------------------------------------------------------------------
def normalize_covariance_matrix(covar):
    """
    Return the correlation matrix from a covariance matrix

    Parameters
    ----------
    covar : numpy.ndarray, shape (N, nBins)
        the covariance matrix to normalize
    """
    N, nBins = covar.shape
    corr = covar.copy()
    
    # normalize the covariance matrix now
    variances = np.diag(corr)
    i, j = np.triu_indices(nBins)
    corr[i, j] /= np.sqrt(variances[i])*np.sqrt(variances[j])
    
    i, j = np.tril_indices(nBins, k=-1)
    corr[i, j] /= np.sqrt(variances[i])*np.sqrt(variances[j])
    
    return corr
#end normalize_covariance_matrix

#-------------------------------------------------------------------------------
def compute_covariance_matrix(X):
    """
    Computes the maximum-likelihood covariance matrix of the data given by X.
    Simple wrapper of the numpy.cov function. See "Estimation of covariance 
    matrices" on Wikipedia, or Barlow 1991 for proof.
    
    Parameters
    ----------
    X : numpy.ndarray, shape (N, N_bins)
        the data array where each column represents a variable and each
        row represents the observations over all variables
    
    Returns
    -------
    covar : numpy.ndarray, shape (N_bins, N_bins)
        the covariance matrix
    """
    return np.cov(X, rowvar=0)
#end compute_covariance_matrix

#-------------------------------------------------------------------------------
def mean_correlations(corr_matrix):
    """
    Compute the mean correlations between bins of a given separation, from 
    a correlation matrix
    
    Parameters
    ----------
    corr_matrix : numpy.ndarray, shape (N_bins, N_bins)
        the correlation matrix
        
    Returns
    -------
    mean_corr : numpy.ndarray, shape (N_bins, )
        the mean correlation separated by a given distance  
    """
    N_bins = corr_matrix.shape[0]
    x = collections.defaultdict(list)
    for i in xrange(N_bins):
        for j in xrange(i, N_bins):
            sep = abs(i-j)
            x[sep].append(corr_matrix[i, j])
            
    mean_corr = np.zeros(N_bins)
    for sep, vals in x.iteritems():
        mean_corr[sep] = np.mean(vals)
    return mean_corr
#end mean_correlations

#-------------------------------------------------------------------------------
def compute_delta_chisq(data, covar_matrix, model, return_chi_sqs=False):
    """
    Compute the square root of chisq_null - chisq_model
    
    Parameters
    ----------
    data : numpy.ndarray
        the array holding the data points
    covar_matrix : numpy.ndarray
        the covariance matrix of the data
    model : numpy.ndarray
        the expected model data points
    return_chi_sqs : bool, optional
        if True, also return chisq_null and chisq_model as the 2nd 
        and 3rd arguments; default = False
    """
    C_inv = np.linalg.inv(covar_matrix)
    chisq_null = np.dot(data, np.dot(C_inv, data))
    chisq_model = np.dot(data-model, np.dot(C_inv, data-model))
    
    sig = np.sqrt(abs(chisq_null - chisq_model))
    if return_chi_sqs:
        return sig, chisq_null, chisq_model
    else:
        return sig
#end compute_delta_chisq

#-------------------------------------------------------------------------------
def compute_null_significance(data, covar_matrix, N_trials=1e6):
    """
    Generate random correlated deviates and compute the significance of the 
    data away from a null signal. Assuming normal distribution errors, return
    the p-value and the sigma value
    
    Parameters
    ----------
    data : numpy.ndarray
        the data array
    covar_matrix : numpy.ndarray
        the covariance matrix
    N_trials : float, optional
        the number of trials to do; default = 1e6
        
    Returns
    -------
    p_value : float
        the probability that random variates have this chi-squared
    sigma : float
        the significance in sigma, assuming gaussian statistics
    """
    nBins = len(data)
    N_trials = int(N_trials)
    
    # do the cholesky decomposition, so C =  A*A.T.conj()
    A = np.linalg.cholesky(covar_matrix)
    
    # get the chisq of the data
    C_inv = np.linalg.inv(covar_matrix)
    chi2_data = np.dot(data, np.dot(C_inv, data))
    
    # compute the correlated data
    devs_uncorrelated = np.vstack([np.random.normal(loc=0., size=N_trials) for i in range(nBins)]).T
    correlated_data = np.zeros(devs_uncorrelated.shape)
    for i in xrange(N_trials):
        correlated_data[i,:] = np.dot(A, devs_uncorrelated[i,:])
        
    # count number of trials where chi2 > chi2_data
    N_larger = 0 
    for i in xrange(int(N_trials)):
        chi2 = np.dot(correlated_data[i,:], np.dot(C_inv, correlated_data[i,:]))
	
        if (chi2 > chi2_data): 
            N_larger += 1
	
    print "N_trials=%d n_larger=%d chi2_data=%f" %(N_trials, N_larger, chi2_data)
    
    p_value = 1.0*N_larger/N_trials
    sigma = getSigmaFromPValue(p_value)
    
    return p_value, sigma, chi2_data
#end compute_null_significance

#-------------------------------------------------------------------------------                                                               
def getSigmaFromChiSquared(chi_sq, dof):
    """
    Compute the significance (in sigma) of a chi squared value and degrees 
    of freedom (from Numerical Recipes, section 6.2)
    
    Parameters
    ----------
    chi_sq : float 
        the chi squared value
    dof : int 
        the degrees of freedom
    
    Returns
    -------
    sigma : float
        the significance in sigma
    p-value : float
        p_value is the probability of obtaining a test statistic at least as 
        extreme as the one that was actually observed
    """
    def func1(x):
        return x - 1.0 + gammaincc(dof/2.0, chi_sq/2.0)
    
    # first calculate p, the confidence limit
    p = sopt.bisect(func1, 0, 1.0)

    def func2(x):
        val = sint.quad(lambda y: 1.0/np.sqrt(2.0*np.pi)*np.exp(-y**2/2.0), -x, x)
        return val[0] - p

    # now calculate sigma
    sigma = sopt.bisect(func2, 0, 100)

    # p is probability that random variates could have this chi-squared value
    return sigma, 1-p
#end getSigmaFromChiSquared
    
#-------------------------------------------------------------------------------
def getSigmaFromPValue(p_value):
    """
    Compute the significance (in sigma) of a given p-value
    
    Parameters
    ----------
    p_value : float 
        p_value is the probability of obtaining a test statistic at least as 
        extreme as the one that was actually observed

    Returns
    -------
    sigma : float
        the significance in sigma
    """
    def objective(x):
        return 1. - erf(x/np.sqrt(2.)) - p_value

    # now calculate sigma
    sigma = sopt.bisect(objective, 0, 100)

    return sigma
#end getSigmaFromPValue

#-------------------------------------------------------------------------------
