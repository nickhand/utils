"""
 gaussian_process.py
 class to implement gaussian process regression, as adapted and simplified 
 from GaPP module (http://www.acgc.uct.ac.za/~seikel/GAPP/)
 
 author: Nick Hand
 contact: nhand@berkeley.edu
 creation date: 05/12/2013
"""
import numpy as np
from scipy import optimize, rand
from numpy import linalg
import tools
import covariance, model
import warnings, copy  


class paramlist(list):
    """
    @brief a class that behaves as a list that allows field access using 
           the input attributes. 
    """
    def __init__(self, *args):
        """
        @brief args are the names of the fields, which must be provided
        """
        if len(args) == 0: raise Exception("No field names provided; "
                                           "paramlist will be empty forever")
        self.dict_ = {}
        for key in args: self.dict_[key] = None
        for i in range(len(args)): self.append(None)
        
    @property
    def dict(self):
        return self.dict_
        
    @property
    def shape(self):
        return np.array(self).shape
        
    @property
    def size(self):
        return self.__len__()
        
    def keys(self):
        return self.dict_.keys()
        
    def items(self):
        return self.dict_.items()
        
    def values(self):
        return self.dict_.values()
    
    def iteritems(self):
        return self.dict_.iteritems()
    
    def copy(self):
        return copy.deepcopy(self)
        
    def __getitem__(self, key):
        if type(key) is str:
            return dict.__getitem__(self.dict_, key)
        elif np.isscalar(key):
            return list.__getitem__(self, int(key))
        else:
            raise TypeError("indices must be scalar or str, not %s" %type(key))
    
    def __setitem__(self, key, value):
        
        if type(key) is str:
            try:
                index = np.where(np.array(self.keys()) == key)[0][0]
                list.__setitem__(self, index, value)
            except:
                raise KeyError('%s' %key)
            dict.__setitem__(self.dict_, key, value)
        elif np.isscalar(key):
            list.__setitem__(self, int(key), value)
            dict.__setitem__(self.dict_, self.keys()[key], value)

        else:
            raise TypeError("indices must be scalar or str, not %s" %type(key))
        
    def __str__(self):
        return "{" + ", ".join("%r: %r" 
                        %(key, self.dict_[key]) for key in self.keys()) + "}"
    def __repr__(self):
        return self.__str__()

#-------------------------------------------------------------------------------
class GaussianProcess(model.model):
    """
    @brief the Gaussian Process class

    Parameters
    ----------
    covfunction : string or callable, optional
        A stationary autocorrelation function returning the autocorrelation
        between two points x and x'.
        Default assumes a squared-exponential covariance function.
        Built-in correlation models are::

            'squared_exponential', 'matern32', 'matern52'

    theta0 : double, array_like or dict
        (Initial) values of the hyperparameters of the covariance function. 
        theta[0] denotes the signal variance and theta[1] the length scale, and 
        theta[2] denotes the noise scale. Can also be dict with keys 
        "signal_variance", "length_scale", and "noise_scale f
        "
        
    mu : callable, optional
        A priori mean function that is subtracted from the data before 
        the GP is performed. Additional args passed through muargs. 
    
    prior : callable, optional
        Prior on the hyperparameters theta. First argument must be theta and it
        must return non-negative number. Additional parameters passed through
        priorargs.
        
    gradprior : callable, optional
       Function that returns a tuple containing the gradient of prior with 
       respect to theta. gradprior needs to be provided if prior != None 
       and grad == True. Additional parameters passed through priorargs.
      
    grad : boolean, optional
       True if the gradient of the covariance function is to be used to train 
       the hyperparameters. This training method is faster than the alternative 
       method. If grad == True and prior != None, one needs to provide the 
       gradient of the prior gradprior. Default is True
       
    verbose : boolean, optional
        A boolean specifying the verbose level.
        Default is verbose = False.

    random_start : int, optional
        The number of times the Maximum Likelihood Estimation should be
        performed from a random starting point.
        The first MLE always uses the specified starting point (theta0),
        the next starting points are picked at random according to an
        exponential distribution (log-uniform on [thetaL, thetaU]).
        Default does not use random starting point (random_start = 1).
        
    fixed : boolean, array_like, optional
        Array equal in length to theta0 or len(theta) + 1, indicating which
        hyperparameters should be held fixed. Default is False for 
        'signal_variance' and 'length_scale' and True for 'noise_scale'
        
    bounds : tuple, array_like, optional
        Array of len(theta0) specifying (min, max) for each hyperparameter. 
        Defaults are (0, None) for each hyperparameter
    """
    
    _covariance_types = {
        'squared_exponential': covariance.SquaredExponential,
        'matern32': covariance.Matern32,
        'matern52': covariance.Matern52}
        
    def __init__(self, theta0, covfunction='squared_exponential',
                     mu = None, muargs=(), prior=None, gradprior=None, 
                     priorargs=(), grad=True, verbose=False, random_start=1, 
                     fixed = None, bounds = None):

        self.covf         = covfunction
        self.verbose      = verbose
        self.random_start = random_start
        
        self.mu        = mu
        self.muargs    = muargs   
        self.prior     = prior
        self.gradprior = gradprior
        self.grad      = grad 
        self.priorargs = priorargs

        # initialize the paramlists we need
        pnames = self._get_param_names()
        self.theta0 = paramlist(*pnames)
        self.bounds = paramlist(*pnames)
        self.fixed  = paramlist(*pnames)
        
        self._initialize_paramlist(theta0, self.theta0)
        self._initialize_paramlist(bounds, self.bounds)
        self._initialize_paramlist(fixed, self.fixed)

        # verify the input parameters
        self._check_params()
        
        # initialize the optimal theta values
        self.theta  = self.theta0.copy()
    #end __init__
    
    #---------------------------------------------------------------------------
    def get_signal_variance(self):
        return self.theta[0]
    
    def get_length_scale(self):
        return self.theta[1]
    
    def get_noise_scale(self):
        return self.theta[2]  
    
    signal_variance = property(get_signal_variance)
    length_scale    = property(get_length_scale)
    noise_scale     = property(get_noise_scale)  
    #---------------------------------------------------------------------------
    def _get_param_names(self):
        return ['signal_variance', 'length_scale', 'noise_scale']
    #end _get_param_names
    #---------------------------------------------------------------------------
    def _initialize_paramlist(self, inlist, paramlist):
        if inlist is not None:
            for i, x in enumerate(inlist):
                paramlist[i] = x
    #end _initialize_paramlist
    
    #---------------------------------------------------------------------------
    def _set_data(self, X, Y, Sigma):
        """
        @brief format and set the observational data
        """
        n = len(X)
        # number of data points
        self.n_samples = n
        if (len(Y) == n and len(Sigma) == n):
            
            # format X and Y based on input shapes
            if(np.shape(X) == (n,)):
                X = np.reshape(X, (n, 1))
                
            self.X = np.array(X)
            self.Y = np.array(Y)
            
            # format the input Sigma array based on shape
            if (np.shape(Sigma) == (n, n)):
                self.Sigma = np.array(Sigma)  
            elif (np.shape(Sigma) in [(n,), (n, 1)]):
                # turn vector into diagonal covariance matrix
                self.Sigma = Sigma * np.eye(n) * Sigma
            else:
                raise ValueError("Input Sigma array must be have shape of one "
                                "of the following: [(N, ), (N, 1), (N, N)]")

        else:
            raise ValueError("Input X, Y and Sigma must have same length")
            
        # subtract mean mu prior, if given
        if self.mu is not None:
            self.Y_mu = np.array([self.Y[i] - self.mu(self.X[i], *self.muargs) \
                                    for i in range(n)])
        else:
            self.Y_mu = self.Y[:]
    #end _set_data
    
    #---------------------------------------------------------------------------
    def fit(self, X, Y, Sigma):
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

        Returns
        -------
        gp : self
            A fitted Gaussian Process model object awaiting data to perform
            predictions.
        """

        
        # set the observational data
        self._set_data(X, Y, Sigma)
        
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
            An array with shape (n_eval, ) with the predicted value, f(x)

        sigma : array_like
            An array with shape (n_eval, ) with the standard deviation at x.
        """
        
        n = len(Xstar)
        # Check input shapes
        if(np.shape(Xstar) == (n,)):
            Xstar = np.reshape(Xstar, (n, 1))
            
        n_eval, n_dim_Xstar = Xstar.shape
        n_samples, n_dim_X = self.X.shape

        # Run input checks
        if n_dim_Xstar != n_dim_X:
            raise ValueError(("The number of dimensions in Xstar "
                              "(Xstar.shape[1] = %d) "
                              "should match the sample size used for fit() "
                              "which is %d.") % (n_dim_Xstar, n_dim_X))


        fmean = np.zeros(n)
        fstd = np.zeros(n)
        
        for i in range(n):
            
            thisXstar = Xstar[i, :]
            nstar = thisXstar.shape[0]
            
            # Get pairwise componentwise L1-distances to the input training set
            dx = tools.l1_distances(tools.array2d(thisXstar), self.X)
        
            # the covariance vector between these distances and training set
            kstar = self.covf.covfunc(self.theta, dx).T
            kstar = kstar.flatten()
            
            # the predictive mean
            mean = np.dot(kstar.T, self.alpha)
        
            # calculate predictive standard deviation
            v = linalg.solve(self.L, kstar)
        
            # now compute cov(Xstar, Xstar)
            dxx = tools.l1_distances(tools.array2d(thisXstar))
            covstar = self.covf.covfunc(self.theta, dxx).T
            covstar = covstar.flatten()[0]
        
            var = covstar - np.dot(v.T, v)
            
            if (self.mu != None):
                mean += self.mu(Xstar, *self.muargs)
            
            if var < 0.: var = 0
            
            fmean[i] = mean
            fstd[i] = np.sqrt(var)
            
        return fmean, fstd
    #end predict
    
    #---------------------------------------------------------------------------
    def likelihood_info(self, theta=None):
        """
        @brief this function determines the best linear unbiased prediction
        parameters and evaluates the reduced likelihood function for the given
        autocorrelation parameters theta.

        Maximizing this function wrt the autocorrelation parameters theta is
        equivalent to maximizing the likelihood of the assumed joint Gaussian
        distribution of the observations y evaluated onto the design of
        experiments X.

        Parameters
        ----------
        theta : array_like, optional
            An array containing the autocorrelation parameters at which the
            Gaussian Process model parameters should be determined.
            Default uses the built-in autocorrelation parameters
            (ie theta = self.theta).

        Returns
        -------
        likelihood_function_value : double
            The value of the likelihood function associated to the
            given autocorrelation parameters theta.
            
        grad_likelihood_function_value : double
            The value of the gradient of likelihood function associated to the
            given autocorrelation parameters theta.

        par : dict
            A dictionary containing the requested Gaussian Process model
            parameters:
        """

        if theta is None:
            # Use built-in autocorrelation parameters
            theta = self.theta

        # Initialize output
        likelihood_function_value = - np.inf
        par = {}

        # Retrieve data
        n_samples = self.n_samples
        D = self.D
    
        # Set up covariance function kernel, K
        k = self.covf.covfunc(theta, D)
        
        A = theta[-1]**2 * self.Sigma[:, :] + \
                                    k.reshape(n_samples, n_samples)[:, :]
       
        # calculate alpha = A^{-1} Y_mu using a Cholesky decomposition
        try:
            L = linalg.cholesky(A)
            b = linalg.solve(L, self.Y_mu)
            alpha = linalg.solve(np.transpose(L), b)
        except linalg.LinAlgError: 
            # A not positive definite
            alpha = None
            L = None
            
        par['L'] = L
        par['alpha'] = alpha
        
        # get the likelihood function
        logp = self._likelihood_function(theta, L, alpha)

        # get the derivative of the likelihood function, if we are using it
        if self.grad:
            grad_logp = self._grad_likelihood_function(theta, L, alpha)
        else:
            grad_logp = None
            
        return (logp, grad_logp, par)
    #end likelihood_info
    
    #---------------------------------------------------------------------------
    def _grad_covariance(self, theta):
        """
        @brief internal method for computing the gradient of the covariance 
               matrix with respect to theta
        """
        n = self.n_samples
        gradK = np.zeros((self.nhparams, n, n))
        
        gradcovf = self.covf.gradcovfunc(theta, self.D)
        gradK[:, ] = gradcovf.reshape(-1, n, n)[:, :]
        
        gradscale = np.zeros((1, n, n))
        scale = theta[-1]
        gradscale[0, ] = 2 * scale * np.eye(n) * np.diagonal(self.Sigma)
        gradK = np.concatenate((gradK,gradscale))
            
        return gradK
    #end _grad_covariance
    
    #---------------------------------------------------------------------------
    def _grad_likelihood_function(self, theta, L, alpha):
        """
        @brief internal method to compute the gradient of the log marginal 
               likelihood and its gradient with respect to theta
        """
        
        n = self.n_samples       
        if (alpha == None):
            try:
                self.gradlogp = 0.9 * self.gradlogp
                return (np.array(self.gradlogp))
            except:
                raise ValueError("Invalid hyperparameters."
                                "Covariance matrix not positive definite")
                 
        # gradient of covariance
        gradK = self._grad_covariance(theta)      
                 
        # number of hyperparameters plus 1 for noise scale
        nh = len(gradK)
        
        # calculate trace of alpha*alpha^T gradK
        traaTgradK = np.zeros(nh)
        for t in range(nh):
            aTgradK = np.zeros(n)
            for i in range(n):
                aTgradK[i] = np.sum(alpha[:] * gradK[t, :, i])
            traaTgradK[t] = np.sum(alpha[:] * aTgradK[:])
        
        # calculate trace of A^{-1} gradK
        invL = tools.chol_inv(L)
        invA = np.dot(invL.T, invL)
        trinvAgradK = np.zeros(nh)
        for t in range(nh):
            for i in range(n):
                trinvAgradK[t] = trinvAgradK[t] + np.sum(invA[i, :] * 
                                                            gradK[t, :, i])
        
        # gradient of the prior log likelihood
        gradpriorlogp = np.zeros(nh)
        if (self.gradprior != None):
            gradpriorp = self.gradprior(theta, *self.priorargs)
            if (self.prior == None):
                warnings.warn("No prior given in "
                              "GaussianProcess.grad_likelihood_function. "
                              "gradprior will be ignored")
            else:
                priorp = self.prior(theta, *self.priorargs)
                for t in range(nh):
                    if (priorp == 0.0 and gradpriorp[t] == 0.0):
                        gradpriorlogp[t] = 0.0
                    elif (priorp <= 0.0):
                        gradpriorlogp[t] = np.sign(gradpriorp[t]) * 1.0e20
                    else:
                        gradpriorlogp[t] = gradpriorp[t]/priorp
        
        # gradient of the negative log likelihood
        gradlogp = np.array(-0.5 * (traaTgradK[:] - trinvAgradK[:]) - gradpriorlogp)
        self.gradlogp = gradlogp
        
        return -1.*gradlogp
    #end _grad_likelihood_function
    
    #---------------------------------------------------------------------------
    def _likelihood_function(self, theta, L, alpha):
        """
        @brief internal method to compute the likelihood function 
        log p(y|X,theta), given the hyperparameters theta and matrices L, alpha
        """
        
        if (self.prior != None):
            priorp = self.prior(theta, *self.priorargs)
            if (priorp < 0.0):
                warnings.warn("invalid prior in "
                             "GaussianProcess.likelihood_function. Negative "
                             "prior will be treated as prior = 0")
                return -1.0e20
            if (priorp == 0.0):
                return -1e20
            priorlogp = np.log(priorp)
        else:
            priorlogp = 0.
            
        # compute the log marginal likelihood
        if (alpha == None):
            logp = -(1.0e20 - priorlogp)
        else:
            logp = (-0.5 * np.dot(self.Y_mu.T, alpha) - 
                      np.sum(np.log(np.diagonal(L))) - self.n_samples/2 * \
                      np.log(2*np.pi) + priorlogp)

        return logp
    #end _likelihood_function
    
    #---------------------------------------------------------------------------
    def optimize(self):
        """
        @brief this function optimizes the hyperparameters by maximimizng the 
        the likelihood function. (Minimization of the negative likelihood 
        function is used for convenience)

        Parameters
        ----------
        self : All parameters are stored in the Gaussian Process model object.

        Returns
        -------
        optimal_theta : array_like
            The best set of autocorrelation parameters (the sought maximizer of
            the reduced likelihood function).

        optimal_reduced_likelihood_function_value : double
            The optimal reduced likelihood function value.

        optimal_par : dict
            The best-fit matrices/parameters associated to optimal_theta.
        """
        
        # Initialize output
        best_optimal_theta = []
        best_optimal_rlf_value = []
        best_optimal_par = []

        if self.verbose:
            if self.random_start > 1:
                print str(self.random_start) + " random starts are required."

        percent_completed = 0.
        for k in range(self.random_start):

            if k == 0:
                # Use specified starting point as first guess
                theta0 = self.theta0
            else:
                # Generate a random starting point log10-uniformly
                # distributed between bounds
                thetaL = np.array([b[0] if b[0] is not None else 0. \
                                    for b in self.bounds])
                inds = np.where(thetaL == 0.)[0]
                for i in inds: thetaL[i] = self.theta0[i]/1e15
                thetaU = np.array([b[1] if b[1] is not None else 1e20 \
                                    for b in self.bounds])
                log10theta0 = np.log10(thetaL) \
                    + rand(self.theta0.size).reshape(self.theta0.shape) \
                    * np.log10(thetaU / thetaL)
                theta0 = 10. ** log10theta0
                
                # reset any values that are fixed
                for i, fix in enumerate(self.fixed):
                    if fix: theta0[i] = self.theta0[i]
                
            # all hyperparameters will be trained
            if np.all([not x for x in self.fixed]):
                
                def logpfunc(theta):
                    logp, gradlogp, par = self.likelihood_info(theta)
                    if gradlogp is None:
                        return -logp
                    else:
                        return (-logp, -gradlogp)
                
                # whether to use gradient information or not
                if self.grad:
                    optimal_theta = optimize.fmin_tnc(logpfunc, theta0, 
                                                     bounds = self.bounds,
                                                     messages=self.verbose*8)[0]
                else:
                    constraints = []
                    for i in range(self.bounds.size):
                        constraints.append(lambda th: th[i]-self.bounds[i, 0])
                        constraints.append(lambda th: self.bounds[i, 1]-th[i])
                             
                    optimal_theta = optimize.fmin_cobyla(logpfunc, theta0, 
                                                         constraints)        
                    
            # some hyperparameters will be trained
            else:                              
                # indices of the hyperparameters that are to be trained
                indices = np.where([not x for x in self.fixed])[0]

                # array of the initial values of these hyperparameters
                inith = np.take(np.array(theta0), indices)
                theta = theta0   # initialize theta
                def logpfunc(th):
                    for i in range(len(indices)):
                        theta[indices[i]] = th[i]
                    
                    logp, gradlogp, par = self.likelihood_info(theta)
                    
                    if gradlogp is None:
                        return -logp
                    else:
                        return (-logp, -gradlogp[indices])
                        
                # whether to use gradient information
                if self.grad:
                    b = np.take(self.bounds, indices, axis=0)
                    th = optimize.fmin_tnc(logpfunc, inith, bounds=b, 
                                            messages=8*self.verbose)[0]
                else:
                    constraints = []
                    for i in indices:
                        constraints.append(lambda th: th[i]-self.bounds[i, 0])
                        constraints.append(lambda th: self.bounds[i, 1]-th[i])
                        
                    th = optimize.fmin_cobyla(logpfunc, theta0, constraints)
                
                for i in range(len(indices)):
                    theta[indices[i]] = th[i]
                optimal_theta = theta
            
            optimal_lf_value, optimal_glf_value, optimal_par = \
                                            self.likelihood_info(optimal_theta)
            
            
            # Compare the new optimizer to the best previous one
            if k > 0:
                if optimal_lf_value > best_optimal_lf_value:
                    best_optimal_lf_value = optimal_lf_value
                    best_optimal_par = optimal_par
                    best_optimal_theta = optimal_theta
            else:
                best_optimal_lf_value = optimal_lf_value
                best_optimal_par = optimal_par
                best_optimal_theta = optimal_theta
                
            if self.verbose and self.random_start > 1:
                if (20 * k) / self.random_start > percent_completed:
                    percent_completed = (20 * k) / self.random_start
                    print "%s completed" % (5 * percent_completed)

        if self.verbose:
            print ("")
            print ("Optimized results:")
            print ("   theta = " + str(best_optimal_theta[:self.nhparams]))
            print ("   scale = " + str(best_optimal_theta[-1]))
            print ("   likelihood = %.4e" %best_optimal_lf_value)
            

        return best_optimal_theta, best_optimal_lf_value, best_optimal_par  
    #end optimize
    
    #---------------------------------------------------------------------------
    def _check_params(self):
        """
        @brief verify all the input parameters and set defaults
        """
        
        # check covariance function
        if not callable(self.covf):
            if self.covf in self._covariance_types:
                self.covf = self._covariance_types[self.covf]
            else:
                raise ValueError(("covfunction should be one of %s or " 
                               + "callable, %s was given.")
                               % (self._covariance_types.keys(), self.covf))
        
        # initialize the covariance function
        self.covf = self.covf()
        
        # number of hyperparameters (without noise scale)
        self.nhparams = self.covf.nparams

        # check the initial hyperparameter values
        if len(self.theta0) != self.nhparams + 1: 
            raise ValueError("Covariance function requires %d input theta "
                             "parameters, %d provided" 
                             % (self.nhparams, len(self.theta0)) )
        
        # force grad to be type bool
        self.grad = bool(self.grad)
        
        # check that gradprior is provided if needed
        if (self.grad and self.prior != None) and self.gradprior is None:
            raise ValueError("If grad = True and prior is not None, then "
                             "gradprior must not be None")
            
        # force verbose type to bool
        self.verbose = bool(self.verbose)
        
        # force random_start type to int
        self.random_start = int(self.random_start)
        
        # set the defaults in the fixed and bounds array
        for pname in self._get_param_names():
            if self.fixed[pname] is None:
                if pname is 'noise_scale':
                    self.fixed[pname]  = True
                    if self.theta0[pname] is None: self.theta0[pname] = 1.
                else:
                    self.fixed[pname] = False
            
            # constrain positive
            self.constrain_positive(pname)                                     
    #end _check_params
#-------------------------------------------------------------------------------

