"""
 tools.py
 utility tools to be used by the gaussian_process module
 
 author: Nick Hand
 contact: nhand@berkeley.edu
 creation date: 05/12/2013
"""
import numpy as np
import scipy.linalg
if hasattr(scipy.linalg.lapack, 'dtrtri'):
    dtrtri = scipy.linalg.lapack.dtrtri
else:
    dtrtri = scipy.linalg.lapack.flapack.dtrtri

def l1_distances(X, Y=None):
    """
    @brief computes the componentwise L1 distances between the 
    vectors in X (sum of absolute value of difference in coordinates)
    
    Parameters
    ----------

    X: array_like
        An array with shape (n_samples_Y, n_features_X)
    
    Y : array_like, optional
        An arra with shape (n_samples_Y, n_features_Y)

    Returns
    -------

    dx: array with shape (n_samples_X * n_samples_Y , n_features_X)
        The array of componentwise L1 distances.
        
    """
    if Y is None:
        Y = X.copy()
    
    n_samples_X, n_features_X = X.shape
    n_samples_Y, n_features_Y = Y.shape
    if n_features_X != n_features_Y:
        raise Exception("X and Y should have the same number of features!")
    D = np.abs(X[:, np.newaxis, :] - Y[np.newaxis, :, :])
    dx = D.reshape((n_samples_X * n_samples_Y, n_features_X))
    
    return dx
#end l1_distances

#-------------------------------------------------------------------------------
def array2d(X, dtype=None, order=None):
    """
    @brief returns at least 2-d array with data from X
    """
    return np.asarray(np.atleast_2d(X), dtype=dtype, order=order)
#end array2d

#-------------------------------------------------------------------------------
def chol_inv(L):
    """
    Inverts a Cholesky lower triangular matrix using scipy lapack

    :param L: lower triangular matrix
    :rtype: inverse of L

    """
    return dtrtri(L, lower = True)[0]