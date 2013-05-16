"""
model.py
base class for the gaussian_process class

Created by Nick Hand on 2013-05-14.
Copyright (c) 2013 UC Berkeley. All rights reserved.
"""
import numpy as np
import copy

class model(object):
    
    def copy(self):
        """
        Returns a (deep) copy of the current model
        """
        return copy.deepcopy(self)
    #end copy
    #---------------------------------------------------------------------------
    def constrain_positive(self, which=None):
        """
        @brief set positive constraints.

        Parameters
        ---------
        which -- scalar int/str or list of strings or ints, which 
                 identifies the parameters. If which=None, all parameters will 
                 be constrained
        """
        if which is None:
            for i in self._get_param_names():
                if self.theta0[i] is not None:
                    self.bounds[i] = (self.theta0[i]/1e15, None)
                else:
                    self.bounds[i] = (1e-20, None)
        else:
            if np.isscalar(which):
                if self.theta0[which] is not None:
                    self.bounds[which] = (self.theta0[which]/1e15, None)
                else:
                    self.bounds[which] = (1e-20, None)
            else:
                for i in which:
                    if self.theta0[i] is not None:
                        self.bounds[i] = (self.theta0[i]/1e15, None)
                    else:
                        self.bounds[i] = (1e-20, None)
        #end constrain_positive
        
    #---------------------------------------------------------------------------
    def unconstrain(self, which=None):
        """
        @brief unconstrain the parameters
        """
        if which is None:
            for i in self._get_param_names():
                self.bounds[i] = (None, None)
        else:
            if np.isscalar(which):
                self.bounds[which] = (None, None)
            else:
                for i in which:
                    self.bounds[i] = (None, None)
        #end unconstrain
        
    #---------------------------------------------------------------------------
    def constrain_bounded(self, which, bounds):
        """
        @brief set bounded constraints.

        Parameters
        ---------
        which --  scalar int/str or list of strings or ints
        bounds -- tuple or list of tuples giving the constraints
        """
        if np.isscalar(which):
            self.bounds[which] = bounds
        else:
            for i, pname in enumerate(which):
                self.bounds[pname] = bounds[i]
        #end constrain_bounded
        
    #----------------------------------------------------------------------------
    def constrain_fixed(self, which, value = None):
        """
        @brief constrain the parameters to be fixed
        """
        if np.isscalar(which):
            self.fixed[which] = True
            if value is not None:
                self.theta0[which] = value
                self.theta[which] = value
        else:
            for i, pname in enumerate(which):
                self.fixed[pname] = True
                self.theta0[pname] = value[i]
                self.theta[pname] = value[i]
        #end constrain_fixed
        
    #---------------------------------------------------------------------------
    def __str__(self):
        """
        @brief return a string describing the parameter names and their 
               ties and constraints
        """
        names = self._get_param_names()
        N = len(names)

        header = ['Name','Value','Constraints']
        values = self.theta 
        
        # sort out the constraints
        constraints = ['']*len(names)
        for i, bound in enumerate(self.bounds):
            constraints[i] = '('+str(bound[0])+', '+str(bound[1])+')'
        for i, fix in enumerate(self.fixed):
            if fix: constraints[i] = 'Fixed'
      
        
        # combine and make the string
        values = ['%.5g' % float(v) for v in values]
        max_names = max([len(names[i]) for i in 
                                        range(len(names))] + [len(header[0])])
        max_values = max([len(values[i]) for i in 
                                        range(len(values))] + [len(header[1])])
        max_constraint = max([len(constraints[i]) for i in 
                                    range(len(constraints))] + [len(header[2])])
        
        cols = np.array([max_names, max_values, max_constraint]) + 4
        columns = cols.sum()

        header_string = ["{h:^{col}}".format(h = header[i], col = cols[i]) 
                                                    for i in range(len(cols))]
        header_string = map(lambda x: '|'.join(x), [header_string])
        separator = '-'*len(header_string[0])
        param_string = ["{n:^{c0}}|{v:^{c1}}|{c:^{c2}}".format(n = names[i], 
                                    v = values[i], c = constraints[i], 
                                    c0 = cols[0], c1 = cols[1], 
                                    c2 = cols[2]) for i in range(len(values))]

        # also print log likelihood if we have fit to data already
        if hasattr(self, 'n_samples'):
            log_like, glog_like, par = self.likelihood_info()
            log_like_str = 'Likelihood: %.3e\n' %(log_like)
            glog_like_str = 'Gradients of likelihood: %s\n' %(glog_like)
        else:
            log_like_str = glog_like_str = ''
        
        return log_like_str + glog_like_str + \
                ('\n'.join([header_string[0], separator]+param_string)) + '\n'

