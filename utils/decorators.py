"""
 decorators.py
 useful function decorators
 
 author: Nick Hand
 contact: nhand@berkeley.edu
 creation date: 05/16/2013
"""
import numpy
from contextlib import contextmanager

@contextmanager
def ignored(*exceptions):
    """
    @brief return a context manager that ignores the specified expections if they
           occur in the body of a with-statement.
    
    For example::
    from contextlib import ignored
    
    with ignored(OSError):
        os.remove('somefile.tmp')
    
    This code is equivalent to:
        try:
            os.remove('somefile.tmp')
        except OSError:
            pass
    
    This will be in python 3.4
    """
    try:
        yield
    except exceptions:
        pass
        
#-------------------------------------------------------------------------------
def call_item_by_item(func):
    """
    @brief Decorator for a function such that an array passed to
    it will be executed item by item.  Return value is the
    same type as the input value (list,ndarray,matrix,etc).

    also up-casts integers.
    """
    def new_func(self,val,*args,**kwargs):
        if type(val) in (int,long):
            val = float(val)
        v_array = numpy.asarray(val)
        v_raveled = v_array.ravel()
        retval = numpy.array([func(self,v,*args,**kwargs) for v in v_raveled],
                             dtype = v_array.dtype)
        retval.resize(v_array.shape)
        if type(val)==numpy.ndarray:
            return retval
        else:
            return type(val)(retval)
    return new_func

#-------------------------------------------------------------------------------
def call_as_array(func):
    """
    @brief decorator for a function such that an array passed to
    it will be executed in one step.  Return value is the
    same type as the input value (float,list,ndarray,matrix,etc).
    
    also up-casts integers.
    """
    def new_func(self,val,*args,**kwargs):
        if type(val) in (int,long):
            val = float(val)
        v_array = numpy.asarray(val)
        v_raveled = v_array.ravel()
        retval = func(self,v_raveled,*args,**kwargs)
        numpy.asarray(retval).resize(v_array.shape)
        if type(val)==numpy.ndarray:
            return retval
        else:
            return type(val)(retval)
    return new_func

