"""
 mat52.py
 a Matern kernel with nu = 5/2
 
 author: Nick Hand
 contact: nhand@berkeley.edu
 creation date: 05/13/2013
"""
import numpy as np

class Matern52(object):
    """
    @brief a class to represent a Matern kernel with nu = 3/2
    """
    def __init__(self):
       self.nparams = 2
    #end __init__
    
    #---------------------------------------------------------------------------
    def covfunc(self, theta, d):
        r"""
        @brief the Matern kernel with \nu = 5/2, given by


        k(x, x_\star) = \sigma_f^2 \mathrm{exp} \left
                        [ -\frac{\sqrt{7}|x - x_\star|}{\ell} \right ] 
                            \left ( 1 +  \frac{\sqrt{7}|x - x_\star|}{\ell} + 
                            \frac{14(x - x_\star)^2}{5\ell^2} + 
                            \frac{7\sqrt{7}|x - x_\star|^3}{15\ell^3} \right )
           
           
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
                
        r = np.sqrt(np.sum(d**2, axis=1))
        covariance = sigmaf**2 * (1 + np.sqrt(5.) * r/l + 5./3. * (r/l)**2) * \
            np.exp(-np.sqrt(5.) * r/l)
        return covariance
    #end covfunc
    
    #---------------------------------------------------------------------------
    def gradcovfunc(self, theta, d):
        """
        @brief the gradient of the squared exponential with respect to the 
        hyperparameters (d/dsigmaf, d/dl)k
        """
        sigmaf, l = theta[:2]
        
        r = np.sqrt(np.sum(d**2, axis=1))
        dk_dsigmaf = 2 * sigmaf * (1 + np.sqrt(5.) * r/l + 5./3. * \
                                         (r/l)**2) * np.exp(-np.sqrt(5.) * r/l)
        dk_dl = 5/3. * sigmaf**2/l**3 * (1 + np.sqrt(5.) * r/l) * r**2 * \
                      np.exp(-np.sqrt(5.) * r/l)
        grad = np.array([dk_dsigmaf, dk_dl])
        return grad
    #end gradcovfunc

    