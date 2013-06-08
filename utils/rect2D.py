# 
#  rect2D.py
#  <for use in a kd-tree module>
#  
#  Created by Nick Hand on 2012-11-25.
# 
import numpy

class point2D(object):
    
    def __init__(self, x, y):
        
        self.xcoord = x
        self.ycoord = y
        
        return 
    
    def x(self):
        """
        @brief the x coordinate
        """
        
        return self.xcoord
        
    def y(self):
        """
        @brief the y coordinate 
        """
        
        return self.ycoord
     
    
class rect2D(object):
    """
    @breif a class that represents an axis-aligned rectangle in the plane
    """
    
    def __init__(self, xmin, ymin, xmax, ymax):
        """
        @brief: construct the rectangle [xmin, xmax] x [ymin, ymax]
                throw an exception if (xmin > xmax) or (ymin > ymax)
        """
              
        if ((xmin >= xmax) or (ymin >= ymax)):
            print xmin, xmax, ymin, ymax
            raise ValueError("Incorrect input values for rectangle coordinates")
        
        self._pLL = point2D(xmin, ymin) # the lower left corner of the rectangle 
        self._pUR = point2D(xmax, ymax) # the upper right corner of the rectangle
        
    def __str__(self):
        """
        @brief: string representation
        """
        return "[%.2f, %.2f] x [%.2f, %.2f]" %(self.xmin(), self.xmax(), self.ymin(), self.ymax())
            
    def xmin(self):
        """
        @brief: minimum x-coordinate of rectangle
        """
        return self._pLL.x()
    
    def ymin(self):
        """
        @brief minimum y-coordinate of rectangle
        """
        return self._pLL.y()
        
    def xmax(self):
        """
        @brief:  maximum x-coordinate of rectangle
        """    
        return self._pUR.x()

    def ymax(self):
        """
        @brief: maximum y-coordinate of rectangle
        """
        return self._pUR.y()
   
    
    def contains(self, p):
        """
        @brief: does this rectangle contain the point p 
                (either inside or on boundary)?
        @param p: point2D object to test on
        """
        xtest = (p.x() >= self.xmin() and p.x() <= self.xmax())
        ytest = (p.y() >= self.ymin() and p.y() <= self.ymax())
        
        return (xtest and ytest)
        
    
    def intersects(self, that):
        """
        @brief: does this rectangle intersect that rectangle 
                (at one or more points)?
        @param that: rect2D object to test on
        """        
        if ((self.xmax() < that.xmin()) or (self.xmin() > that.xmax())):
            return False
        if ((self.ymax() < that.ymin()) or (self.ymin() > that.ymax())):
            return False
        
        return True
    
     
    def distanceTo(self, p):
        """
        @brief: Euclidean distance from point to the closest point in rectangle
        @param p: point2D object to calculate distance to
        """
        # distance = 0 if p is within rectangle
        if (self.contains(p)): 
            return 0. 
        
        d = self.distance2To(p) 
        
        return numpy.sqrt(d)
    
    def distance2To(self, p): 
        """
        @brief: square of Euclidean distance from point p to closest point in rectangle
        @param p: point2D object to calc dist to
        """
    
        # distance = 0 if p is within rectangle
        if (self.contains(p)): 
            return 0. 
        
        dx = 0. 
        dy = 0.
     
        # point lies left of rect
        if (p.x() < self.xmin()): 
            dx = (p.x() - self.xmin())**2
        # point lies right of rect
        elif (p.x() > self.xmax()):
            dx = (p.x() - self.xmax())**2
     
        # point lies below rect
        if (p.y() < self.ymin()):  
            dy = (p.y() - self.ymin())**2
        # point lies above rect
        elif (p.y() > self.ymax()):
            dy = (p.y() - self.ymax())**2
                
        return (dx + dy)
    
    def __eq__(self, that):
        """
        @brief: does this rectangle equal that?
        @param that: rect2D to test equality
        """   
        
        # check null
        if (that is None): 
            return False 
        
        # must be same class
        if not isinstance(that, self.__class__):
            return False
        
        # check corners
        if (self.xmin() != that.xmin()): 
            return False
        if (self.xmax() != that.xmax()): 
            return False
        if (self.ymin() != that.ymin()):
             return False
        if (self.ymax() != that.ymax()):
             return False
        
        return True 
    
    def __hash__(self):
        """
        @brief: the python default hash method
        """
        return hash( (self.xmin(), self.ymin(), self.xmax(), self.ymax()) )
