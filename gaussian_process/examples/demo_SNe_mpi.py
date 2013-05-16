"""
 demo_SNe_mpi.py
 a Gaussian Process regression example using supernova data, fitting in parallel
 using mpi4py
 
 author: Nick Hand
 contact: nhand@berkeley.edu
 creation date: 05/14/2013
"""
import pylab as pl
import numpy as np
import argparse
import sys, os
sys.path.append("..")
from gaussian_process_mpi import GaussianProcessMPI
from glob import glob
from mpi4py import MPI
comm = MPI.COMM_WORLD

def run_demo(args):
    """
    @brief a Gaussian Process regression example that fits several supernovae
           spectra in parallel using mpi4pys
    """
    
    myrank = comm.Get_rank()
    nprocs = comm.Get_size()
    
    # read in the relevant files in correct order
    f1 = glob("../data/SN2011fe/11feM*")
    f1.sort(reverse=True)
    f2 = glob("../data/SN2011fe/11feP*")
    f2.sort()
    files = f1+f2
    
    pl.ion()
    
    for N, f in enumerate(files):
        
        file_root = os.path.splitext(os.path.basename(f))[0]
        time = float(file_root[-3:])/10.
        if 'M' in file_root: time *= -1.
        
        # load the data from inputdata.txt
        X, Y, Yerr = np.loadtxt(f, unpack=True)
  
        n_eval = len(X)
        batch_size = args.batch_size
        resolution = args.resolution
    
        # instanciate a Gaussian Process model, allowing all params to vary
        gp = GaussianProcessMPI(theta0 = [1e-13, 2.0, 1e-13], 
                             covfunction=args.covariance, verbose=True, 
                             fixed=[False, False, False])

        xmin = np.amin(X)
        xmax = np.amax(X)
        nstar = len(X)*resolution
    
        # mesh the input space for evaluations of the prediction
        x = np.linspace(xmin, xmax, nstar)

        # fit to data using Maximum Likelihood Estimation of the parameters
        gp.fit(X, Y, Yerr, batch_size=200)

        # make the prediction on the meshed x-axis
        y_pred, sigma = gp.predict(x)
        
        # plot if you are the master
        if myrank == 0:
        
            # plot the function, the prediction and the 95% confidence interval based on
            # the standard deviation
            pl.cla()
            pl.plot(X, Y, label='Observations')
            pl.plot(x, y_pred, label='Prediction')
            pl.fill(np.concatenate([x, x[::-1]]),
                   np.concatenate([y_pred - 1.9600 * sigma,
                                  (y_pred + 1.9600 * sigma)[::-1]]),
                   alpha=0.5, fc='DarkGoldenRod', ec="None", 
                   label='95% confidence interval')
                
            pl.xlabel(r'$\lambda \ (\AA)$', fontsize=16)
            pl.ylabel('$\mathrm{Flux \ (erg/s/cm^2/\AA)}$', fontsize=16)
            pl.legend(loc='upper right')
            pl.title("SNe 2011fe %+.1f days relative to B-band max" %time)
            pl.ylim(-0.2e-12, 1.2e-12)
            pl.savefig("figures/SN11fe_mpi_%02d.png" %N)
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
    
    if comm.Get_rank() == 0:
        # initialize the output dir
        if not os.path.exists('figures'):
            os.makedirs('figures')
        
    # run the demo
    run_demo(args)
    

