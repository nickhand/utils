"""
 samplinator.py
 class that provides sampling and B-spline fitting of methods
 
 original author: Jake vanderPlas
 editing author: Nick Hand
 contact: nhand@berkeley.edu
 creation date: 05/01/2013
"""

import numpy
from scipy import integrate, interpolate

class with_sampleable_methods:
    """
    Class which allows sampling and B-Spline fitting of class methods
    in derived classes.

    Example:

    #---------------------------------------------------------------
    class foo(with_sampleable_methods):
        def bar(self,x):
            return numpy.sin(x)
        #--
    #--

    F = foo()
    print F(4) #evaluates the function
    F.sample('bar',numpy.arange(0,10,0.1))
    print F(4) #evalueates a pre-computed B-spline of the function
    #---------------------------------------------------------------
    """
    class sampled_function:
        def __init__(self,func,x,*args,**kwargs):
            self.func = func
            self.x = x

            #assign function name
            if func.__name__ != None:
                self.__name__ = self.func.__name__ + \
                                " [Sampled to %i pts]" % len(x)
            else:
                self.__name__ = None

            #assign function doc string
            if func.__doc__ != None:
                self.__doc__ = "Sampled Function : \n\n" + self.func.__doc__
            else:
                self.__doc__ = None

            #set up the b-spline
            try:
                self.tck = interpolate.splrep(x,func(x,*args,**kwargs),s=0)
            except:
                self.tck = interpolate.splrep(x,func,s=0)
        def __call__(self,y):
            return interpolate.splev(y,self.tck,der=0)
    ###
        
    def sample(self,methodname,xrange,*args,**kwargs):
        if not hasattr(self,methodname):
            raise ValueError, methodname

        if self.is_sampled(methodname):
            self.unsample(methodname)

        tmp = getattr(self,methodname)





        setattr(self,methodname,
                with_sampleable_methods.sampled_function(tmp,xrange,
                                                         *args,**kwargs ) )

    def unsample(self,methodname):
        if self.is_sampled(methodname):
            tmp = getattr(self,methodname).func
            setattr( self,methodname,tmp )
        else:
            raise ValueError, "cannot unsample %s" % methodname
        


    def is_sampled(self,methodname):
        return getattr(self,methodname).__class__ == \
               with_sampleable_methods.sampled_function



