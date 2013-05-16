"""
 demo1.py
 a Gaussian Process regression example
 
 author: Nick Hand
 contact: nhand@berkeley.edu
 creation date: 05/14/2013
"""
import pylab as pl
import numpy as np
import argparse
import sys
sys.path.append("..")
from gaussian_process import GaussianProcess

def run_demo(args):
    """
    @brief a Gaussian Process regression example using an input covariance model
    """

    np.random.seed(1)

    def f(x):
        """
        @brief the function to predict.
        """
        return x * np.sin(x)
  
    # the inpute data points
    X = np.linspace(0.1, 9.9, 20)

    # make the observations with added noise
    y = f(X).ravel()
    dy = 0.5 + 1.0 * np.random.random(y.shape)
    noise = np.random.normal(0, dy)
    y += noise
    
    # mesh the input space for evaluations of the prediction
    x = np.linspace(-2, 12, 2*len(X))

    # instanciate a Gaussian Process model, allowing all params to vary
    gp = GaussianProcess(theta0 = [0.5, 2.0, 1.0], 
                         covfunction=args.covariance, verbose=True, 
                         fixed=[False, False, False], random_start=10)

    # fit to data using Maximum Likelihood Estimation of the parameters
    gp.fit(X, y, dy)

    # make the prediction on the meshed x-axis
    y_pred, sigma = gp.predict(x)
    
    # plot the function, the prediction and the 95% confidence interval based on
    # the standard deviation
    fig = pl.figure()
    pl.plot(x, f(x), 'r:', label=r'$f(x) = x \ \mathrm{sin}(x)$')
    pl.errorbar(X.ravel(), y, dy, label='Observations')
    pl.plot(x, y_pred, label='Prediction')
    pl.fill(np.concatenate([x, x[::-1]]),
            np.concatenate([y_pred - 1.9600 * sigma,
                           (y_pred + 1.9600 * sigma)[::-1]]),
            alpha=.2, fc='DarkGoldenRod', ec="None", label='95% confidence interval')
    
    pl.xlabel('$x$', fontsize=16)
    pl.ylabel('$f(x)$', fontsize=16)
    pl.ylim(-15, 20)
    pl.legend(loc='upper left')
    
    pl.show()
#end run_demo

#-------------------------------------------------------------------------------
if __name__ == "__main__":
    
    # parse the command line arguments
    parser = argparse.ArgumentParser(description="example 1 of a Gaussian "
                                     "Process regression")
    parser.add_argument("--covariance", type=str, 
                        choices=['squared_exponential',
                        'matern32', 'matern52'], 
                        default='squared_exponential',
                        help="the covariance model to use")                      
    args = parser.parse_args()
    
    # run the demo
    run_demo(args)
    
    

