"""
 demo_SNe.py
 a Gaussian Process regression example using supernova data
 
 author: Nick Hand
 contact: nhand@berkeley.edu
 creation date: 05/14/2013
"""
import pylab as pl
import numpy as np
import argparse
import sys, os
sys.path.append("..")
from gaussian_process import GaussianProcess
from glob import glob

def run_demo(args):
    """
    @brief a Gaussian Process regression example that fits several supernovae
           spectra
    """
    
    # read in the relevant files in correct order
    f1 = glob("../data/SN2011fe/11feM*")
    f1.sort(reverse=True)
    f2 = glob("../data/SN2011fe/11feP*")
    f2.sort()
    files = f1+f2
    
    pl.ion()
    
    # fit each supernova spectrum in serial
    for i, f in enumerate(files):
        
        file_root = os.path.splitext(os.path.basename(f))[0]
        time = float(file_root[-3:])/10.
        if 'M' in file_root: time *= -1.
        
        # load the data from inputdata.txt
        X, Y, Yerr = np.loadtxt(f, unpack=True)
  
        # save these for later
        X_tot    = X.copy()
        Y_tot    = Y.copy()
        Yerr_tot = Yerr.copy()
    
        n_eval = len(X)
        batch_size = args.batch_size
        resolution = args.resolution
    
        xrec_full = []
        yrec_full = []
        yerr_full = []
    
        # instanciate a Gaussian Process model, allowing all params to vary
        gp = GaussianProcess(theta0 = [1e-13, 2.0, 1e-13], 
                             covfunction=args.covariance, verbose=True, 
                             fixed=[False, False, False])
                         
        nbatches = max(1, n_eval / batch_size + 1)
        # fit the spectra in batches along the x axis
        for k in range(nbatches):

            batch_from = k * batch_size 
            batch_to = min([(k + 1) * batch_size + 1, n_eval + 1])
            if k == nbatches-1: batch_to = len(X_tot)
            
            xmin = np.amin(X_tot[batch_from:batch_to])
            xmax = np.amax(X_tot[batch_from:batch_to])
            nstar = len(X_tot[batch_from:batch_to])*resolution
        
            batch_to += 0.1*batch_size
            batch_from -= 0.1*batch_size
            if batch_from < 0: batch_from = 0
        
            X  = X_tot[batch_from:batch_to]
            Y  = Y_tot[batch_from:batch_to]
            Yerr = Yerr_tot[batch_from:batch_to]
    
            # mesh the input space for evaluations of the prediction
            x = np.linspace(xmin, xmax, nstar)

            # fit to data using Maximum Likelihood Estimation of the parameters
            gp.fit(X, Y, Yerr)

            # make the prediction on the meshed x-axis
            y_pred, sigma = gp.predict(x)
    
    
            xrec_full += list(x)
            yrec_full += list(y_pred)
            yerr_full += list(sigma)
        
        yerr_full = np.array(yerr_full)
        yrec_full = np.array(yrec_full)
    
        # plot the function, the prediction and the 95% confidence interval based on
        # the standard deviation
        pl.cla()
        pl.plot(X_tot, Y_tot, label='Observations')
        pl.plot(xrec_full, yrec_full, label='Prediction')
        pl.fill(np.concatenate([xrec_full, xrec_full[::-1]]),
               np.concatenate([yrec_full - 1.9600 * yerr_full,
                              (yrec_full + 1.9600 * yerr_full)[::-1]]),
               alpha=0.5, fc='DarkGoldenRod', ec="None", 
               label='95% confidence interval')
    
        pl.xlabel(r'$\lambda \ (\AA)$', fontsize=16)
        pl.ylabel('$\mathrm{Flux \ (erg/s/cm^2/\AA)}$', fontsize=16)
        pl.legend(loc='upper right')
        pl.title("SNe 2011fe %+.1f days relative to B-band max" %time)
        pl.ylim(-0.2e-12, 1.2e-12)
        pl.savefig("figures/SN11fe_%02d.png" %i)
        pl.draw()
#end run_demo

#-------------------------------------------------------------------------------
if __name__ == "__main__":
    
    # read in the command line arguments
    parser = argparse.ArgumentParser(description="example 1 of a Gaussian "
                                                 "Process regression")
    parser.add_argument("--covariance", type=str, 
                        choices=['squared_exponential',
                        'matern32', 'matern52'], 
                        default='squared_exponential',
                        help="the covariance model to use")
    parser.add_argument("--batch_size", type=int, 
                        default=200,
                        help="the x-width to fit spectra at one time")
    parser.add_argument("--resolution", type=int, 
                        default=2,
                        help="the factor increase in number of pts on x-axis")
    args = parser.parse_args()
    
    # initialize the output dir
    if not os.path.exists('figures'):
        os.makedirs('figures')
        
    # run the demo
    run_demo(args)
    

