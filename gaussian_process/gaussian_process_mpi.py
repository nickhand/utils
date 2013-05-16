"""
 gaussian_process_mpi.py
 class to implement gaussian process regression in parallel using mpi4py
 
 author: Nick Hand
 contact: nhand@berkeley.edu
 creation date: 05/12/2013
"""
import numpy as np
from numpy import linalg
import tools
import copy
from gaussian_process import GaussianProcess
from mpi4py import MPI
comm = MPI.COMM_WORLD

class GaussianProcessMPI(GaussianProcess):
    
    def fit(self, X, Y, Sigma, batch_size=200, overlap=0.2):
        """
        @brief the method whether the fitting magic happens

        Parameters
        ----------
        X : double array_like
            An array with shape (n_samples, ) or (n_samples, n_dim) with the 
            input at which observations were made.

        Y : double array_like
            An array containing the observations f(X)
        
        Sigma : double array_like
            An array of shape (n_samples, ) containing the measurement errors 
            of Y, or a (n_samples, n_samples) covariance matrix of the data.
            
        batch_size : double, optional
            the width in points along the x-axis to fit one at a time
            
        overlap : double, optional
            the percentage of points to overlap between batches that are fit

        Returns
        -------
        gp : self
            A fitted Gaussian Process model object awaiting data to perform
            predictions.
        """
        nprocs = comm.Get_size()
        myrank = comm.Get_rank()
        
        nbatches = max(1, len(X)/batch_size)
        self.nbatches = nbatches
        
        n_samples = len(X)
        self.fits = []
        self.xbounds = []
        
        # do given batch of the fit based on rank
        for k in range(myrank, nbatches, nprocs):
            
            batch_from = k * batch_size 
            batch_to = min([(k + 1) * batch_size + 1, n_samples + 1])
            if k == nbatches-1: batch_to = len(X)
        
            # keep track of the x bounds
            xmin = np.amin(X[batch_from:batch_to])
            xmax = np.amax(X[batch_from:batch_to])
            self.xbounds.append((xmin, xmax))
            
            batch_to += overlap*batch_size
            batch_from -= overlap*batch_size
            if batch_from < 0: batch_from = 0
        
            X_cut    = X[batch_from:batch_to]
            Y_cut    = Y[batch_from:batch_to]
            Sigma_cut = Sigma[batch_from:batch_to]
            
            # set the observational data
            self._set_data(X_cut, Y_cut, Sigma_cut)
        
            # calculate matrix of distances D between input X samples
            D = tools.l1_distances(self.X)
            self.D = D

            # do the fitting
            # Maximum Likelihood Estimation of the parameters
            if self.verbose:
                print("Performing Maximum Likelihood Estimation of the "
                    + "autocorrelation parameters...")
                
            self.theta, self.likelihood_function_value, par = self.optimize()
            if np.isinf(self.likelihood_function_value):
                raise Exception("Bad parameter region. Try increasing upper bound")
        
            self.alpha = par['alpha']
            self.L     = par['L']
            
            self.fits.append(copy.deepcopy(self))
        
        return self
    #end fit
    #---------------------------------------------------------------------------
    def predict(self, Xstar):
        """
        @brief evaluate the Gaussian Process model at X.

        Parameters
        ----------
        Xstar : array_like
            An array with shape (n_eval, n_features) giving the point(s) at
            which the prediction(s) should be made.

        Returns
        -------
        y : array_like
            An array with shape (n_eval, ) with the Best Linear Unbiased
            Prediction at x.

        sigma : array_like
            An array with shape (n_eval, ) with the standard deviation at x.
        """
        nprocs = comm.Get_size()
        myrank = comm.Get_rank()
    
        data_pred = [] # the predicted data
        
        # do given batch of the fit based on rank
        cnt = 0
        for k in range(myrank, self.nbatches, nprocs):
            
            # trim Xstar
            inds = np.where((Xstar>=self.xbounds[cnt][0])*(Xstar<=self.xbounds[cnt][1]))
            Xstar_cut = Xstar[inds]
            
            # check for Xstar outside fitting X range
            if k == 0:
                inds = np.where(Xstar < self.xbounds[cnt][0])
                Xstar_cut = np.append(Xstar[inds], Xstar_cut)
            if k == self.nbatches-1:
                inds = np.where(Xstar > self.xbounds[cnt][1])
                Xstar_cut = np.append(Xstar_cut, Xstar[inds])
            
            n = len(Xstar_cut)
            
            # Check input shapes
            if(np.shape(Xstar_cut) == (n,)):
                Xstar_cut = np.reshape(Xstar_cut, (n, 1))
            
            n_eval, n_dim_Xstar = Xstar_cut.shape
            n_samples, n_dim_X = self.X.shape

            # Run input checks
            if n_dim_Xstar != n_dim_X:
                raise ValueError(("The number of dimensions in Xstar "
                                  "(Xstar.shape[1] = %d) "
                                  "should match the sample size used for fit() "
                                  "which is %d.") % (n_dim_Xstar, n_dim_X))


            # get the right fit
            fit = self.fits[cnt]
            
            fmean = np.zeros(n)
            fstd = np.zeros(n)
        
            for i in range(n):
            
                thisXstar = Xstar_cut[i, :]
                nstar = thisXstar.shape[0]
            
                # Get pairwise componentwise L1-distances to the input training set
                dx = tools.l1_distances(tools.array2d(thisXstar), fit.X)
        
                # the covariance vector between these distances and training set
                kstar = self.covf.covfunc(fit.theta, dx).T
                kstar = kstar.flatten()
            
                # the predictive mean
                mean = np.dot(kstar.T, fit.alpha)
        
                # calculate predictive standard deviation
                v = linalg.solve(fit.L, kstar)
        
                # now compute cov(Xstar, Xstar)
                dxx = tools.l1_distances(tools.array2d(thisXstar))
                covstar = self.covf.covfunc(fit.theta, dxx).T
                covstar = covstar.flatten()[0]
        
                var = covstar - np.dot(v.T, v)
            
                if (self.mu != None):
                    mean += self.mu(Xstar, *self.muargs)
            
                if var < 0.: var = 0
            
                fmean[i] = mean
                fstd[i] = np.sqrt(var)
            
            cnt += 1
           
            # save the predictions
            data_pred.append((k, Xstar_cut, fmean, fstd))
         
        # gather the reconstructed data to the master   
        data_pred = comm.gather(data_pred, root=0) 
        
        # the master combines and returns
        if myrank == 0: 
            x_rec, y_rec, yerr_rec = self._combine_data(data_pred)
            
            return y_rec, yerr_rec
        else:
            return None, None
    #end predict
    
    #---------------------------------------------------------------------------
    def _combine_data(self, data_pred):
        """
        @brief internal method to recombine the data in data_pred to 
        x, y and sigma arrays  
        """
        
        data_sorted = [None]*self.nbatches
        cnt = 0
        for i in range(comm.Get_size()):
            d = data_pred[i]
            for j in range(len(d)):
                data_sorted[cnt] = d[j]
                cnt += 1

        xrec_full = []
        yrec_full = []
        yerr_full = []
        
        # put the predicted data in the right order
        data_sorted.sort(key=lambda x: x[0])
        
        for (i, x, y_pred, sigma) in data_sorted:
            xrec_full += list(x)
            yrec_full += list(y_pred)
            yerr_full += list(sigma)
            
        return np.array(xrec_full), np.array(yrec_full), np.array(yerr_full)
    #end  _combine_data
    
