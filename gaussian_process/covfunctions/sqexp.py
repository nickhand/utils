"""
 sqexp.py
 a squared exponential covariance function
 
 author: Nick Hand
 contact: nhand@berkeley.edu
 creation date: 05/13/2013
"""
import numpy as np

class SquaredExponential(object):
    """
    @brief a class to represent a squared exponential kernel
    """
    
    # initialize class with initial hyperparameter theta
    def __init__(self):
        
        self.nparams = 2
    #end __init__
    
    #---------------------------------------------------------------------------
    def covfunc(self, theta, d):
        r"""
        Squared exponential covariance function given by

                                                   
            k(x, x_\star) = \sigma_f^2 \mathrm{exp} \
                            \left ( - \frac{(x-x_\star)^2}{2\ell^2} \right )
                                                  

        Parameters
        ----------
        theta : array_like
            An array with length of at least 2 such that \sigma_f = theta[0]
            and \ell = theta[1]

        d : array_like
            An array with shape (n_eval, n_features) giving the componentwise
            distances between locations x and x' at which the covariance model
            should be evaluated.

        Returns
        -------
        r : array_like
            An array with shape (n_eval, ) containing the values of the
            covariance model.
        """    
        sigmaf, l = theta[:2]
        xxl = np.sum((d/l)**2, axis=1)
        covariance = sigmaf**2 * np.exp(-xxl/2.)
        return covariance
    #end covfunc
    
    #---------------------------------------------------------------------------
    def gradcovfunc(self, theta, d):
        """
        @brief the gradient of the squared exponential with respect to the 
        hyperparameters (d/dsigmaf, d/dl)k
        """
        sigmaf, l = theta[:2]   
        xxl = np.sum((d/l)**2, axis=1)
        dk_dsigmaf = 2 * sigmaf * np.exp(-xxl/2.)
        dk_dl = sigmaf**2/l * xxl * np.exp(-xxl/2.)
        grad = np.array([dk_dsigmaf, dk_dl])
        return grad
    #end gradcovfunc
