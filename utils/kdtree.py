"""
 kdtree.py
 kd-tree module to do fast nearest neighbor searches on objects
 
 author: Nick Hand
 contact: nhand@berkeley.edu
 creation date: 06/08/2013
"""
from xcorrSZ import rect2D 
import copy
import numpy as np

class kdtree(object):

    class _Node(object):
        """
        @brief class representing a node of the kd-tree
        """ 
        
        def __init__(self, r, element):  
            """
            @brief: construct a new node with a given point and rectangle
            """    
            
            # the axis-aligned rectangle corresponding to the node
            self._rect = r  
        
            # the left/bottom subtree
            self._lb = None 
            
            # the right/top subtree
            self._rt = None 
        
            # the element being stored in this Node
            self._element = element 
                
            
            return


    def __init__(self, pointClass, xmin, ymin, xmax, ymax):
        """
        @brief: construct an empty kdtree, with total bounding box
        ranging from xmin to xmax and ymin to ymax
        
        @param xmin: min x coordinate (float)
        @param xmax: max x coordinate (float)
        @param ymin: min y coordinate (float)
        @param ymax: max y coordinate (float)
        """
            
        # initially the root Node is None and N = 0
        self._root = None
        self._N = 0
        
        # the bounding box is RA: [0, 360], DEC: [-180, 180]
        self._boundingBox = [xmin, ymin, xmax, ymax]
        self.xmin, self.ymin, self.xmax, self.ymax = xmin, ymin, xmax, ymax
        
        # the list of objects in the tree
        self.objects = []
        
        # save the point class
        self.pointClass = pointClass
        
        return

    def isEmpty(self):
        """
        @brief: return True if the set is empty
        """
        return self._N == 0
        
    def size(self):
        """
        @brief: return the size of the kd-tree
        """
        return self._N
            
    def _insert(self, n, r, p, depth):
        """
        @brief: internal recursive insert method used to insert new node somewhere in 
                the subtree of Node n
        @param n: Node object that is head of tree
        @param r: rect2D object corresponding to the bounding box at curr depth
        @param p: the point to be inserted (point)
        @param depth: current depth in the subtree
            
        @return: the inserted node
        """
            
        # insert node if end of tree is reached
        if (n is None):
            self._N += 1
            self.objects.append(p)
            return kdtreeGeneric._Node(r, p)
        
        # compare Nodes by x coord if on even depth level
        if ((depth % 2) == 0):
            
            # go to the left subtree
            if (p.x() < n._element.x()):
                
                # call insert with left child, keeping track of rectangle bounding box
                if (n._lb is None):
                    rect = rect2D.rect2D(r.xmin(), r.ymin(), n._element.x(), r.ymax())
                    n._lb = self._insert(n._lb, rect, p, depth+1)
                else: 
                    n._lb = self._insert(n._lb, n._lb._rect, p, depth+1)

            # go to the right subtree
            else:
                
                # call insert with right child, keeping track of rectangle bounding box
                if (n._rt is None):
                    rect = rect2D.rect2D(n._element.x(), r.ymin(), r.xmax(), r.ymax())
                    n._rt = self._insert(n._rt, rect, p, depth+1)
                else:
                    n._rt = self._insert(n._rt, n._rt._rect, p, depth+1)
        
        # compare Nodes by y coord if on odd depth level
        else:    
            # go to the left subtree
            if (p.y() < n._element.y()):
                
                # call insert with left child, keeping track of rectangle bounding box
                if (n._lb is None): 
                    rect = rect2D.rect2D(r.xmin(), r.ymin(), r.xmax(), n._element.y())
                    n._lb = self._insert(n._lb, rect, p, depth+1)
                else:
                    n._lb = self._insert(n._lb, n._lb._rect, p, depth+1)
                        
            # go the right subtree
            else:
                
                # all insert with right child, keeping track of rectangle bounding box
                if (n._rt is None):
                    rect = rect2D.rect2D(r.xmin(), n._element.y(), r.xmax(), r.ymax())
                    n._rt = self._insert(n._rt, rect, p, depth+1)
                else:
                    n._rt = self._insert(n._rt, n._rt._rect, p, depth+1)
        
        return n 
        
    def insert(self, obj):
        """
        @brief: add the object specified to the kd-tree
        @param obj: object to be added (sourceSZ.source)
        """
        
        # the initial depth
        depth = 0
        
        # bounding box of root is specified by self.boundingBox
        r = rect2D.rect2D(*self._boundingBox)
        
        # make the point object with necessary functions
        p = self.pointClass(obj.__dict__.copy())
        
        # raise error if ra is not [0, 360] and dec is not [-180, 180]
        if not r.contains(p): raise ValueError('input object %s not in bounding box' %obj)
        
        # insert the object into the tree, recursively
        self._root = self._insert(self._root, r, p, depth)
        
        return
        
    def _contains(self, n , p, depth):
        """
        @brief: internal method used to recursively check if object is contained within the subtree of Node n
        
        @param n: Node object that is head of tree
        @param p: point object to see if is contained
        @param depth: current depth in the subtree
        """

        # if we reach end of tree, return False
        if (n is  None): 
            return False 
        
        # return true if found
        if (n._element == p):
            return True 
        
        # compare Nodes by x coord for even depth level 
        if ((depth % 2) == 0):
            
            # call contains() for left and right child
            if (p.x() < n._element.x()): 
                ans = self._contains(n._lb, p, depth+1)
            else:
                ans = self._contains(n._rt, p, depth+1)
        # compare Nodes by y coord for odd depth level
        else:
            
            # call contains() for left and right child
            if (p.y() < n._element.y()):
                ans = self._contains(n._lb, p, depth+1)
            else:
                ans = self._contains(n._rt, p, depth+1)
        
        # return the answer from the subtree searches
        return ans
    
    
    def contains(self, obj):
        """
        @brief: return True if the tree contains the object specified
        @param obj: the object to test
        
        @return: True/False
        """
        
        depth = 0
        
        # make the point object with necessary functions
        p = self.pointClass(obj.__dict__.copy())
        
        # return the answer of the contains recursive search
        return self._contains(self._root, p, depth) 
    
    def _peek(self, n , p, depth):
        """
        @brief: recursively check if object obj is contained within the subtree of Node n
        @param n: Node object that is head of tree
        @param p: point object to test for
        @param depth: current depth in the subtree
        
        @return: the Node containing the object, if found
        """
        
        # if we reach end of tree, return None
        if (n is  None): 
            return None 
        
        # return the Node containing the object, if found
        if (n._element == p):
            return p 
        
        # compare Nodes by x coord for even depth level 
        if ((depth % 2) == 0):
            
            # call contains() for left and right child
            if (p.x() < n._element.x()): 
                ans = self._peek(n._lb, p, depth+1)
            else:
                ans = self._peek(n._rt, p, depth+1)
        # compare Nodes by y coord for odd depth level
        else:
            
            # call contains() for left and right child
            if (p.y() < n._element.y()):
                ans = self._peek(n._lb, p, depth+1)
            else:
                ans = self._peek(n._rt, p, depth+1)
        
        return ans
    
    def peek(self, obj):
        """
        @brief: return the Node containing obj or None if the tree does not contain obj
        @param obj: object to look for
        
        @return: the copy of obj from the kd-tree or None
        """
        
        depth = 0
        
        # make the point object with necessary functions
        p = self.pointClass(obj.__dict__.copy())
        
        # return the result of the peek search
        return self._peek(self._root, p, depth) 
    
        
    def _rangeSearch(self, n, rect, s, p0, radius):
        """
        @brief: internal recursive method used to find all points in subtree of Node n 
                that fall within the Rect2D object rect
        @param n: Node object that is head of tree to search
        @param rect: rect2D object to search within
        @param s: set object containing all objects found so far 
        @param p0: point2D object to look within
        @param radius: radius to look within
        """
        
        # add current point to set if contained
        if (rect.contains(n._element)):
        
            # make sure its also in the circle then add Node's element
            if n._element.distanceTo(p0) <= radius:
                s.add(n._element)        
        
        if (n._lb is not None):         
            # only explore left subtree if bounding box of child intersects rect
            if (rect.intersects(n._lb._rect)): 
                self._rangeSearch(n._lb, rect, s, p0, radius)           
        
        if (n._rt is not None):
            # only explore left subtree if bounding box of child intersects rect
            if (rect.intersects(n._rt._rect)):     
                self._rangeSearch(n._rt, rect, s, p0, radius)
                
        return
    
    def range(self, obj0, radius, functionToRun=None, checkEdges=True):
        """
        @brief: find all objects in the kd-tree that are inside the circle centered around
                obj0 and within the radius specified
        
        @return s: set object containing all points        
        """        
        
        # make the point object with necessary functions
        p0 = self.pointClass(obj0.__dict__.copy())
        
        # the rectangle to find objects in
        rect = rect2D.rect2D(p0.x()-radius, p0.y()-radius, p0.x()+radius, p0.y()+radius)        
        
        # the set of objects we will return
        s = set()
        
        # if kd-tree is empty, return
        if self.isEmpty(): return s
        
        # check if we are close to edge of bounding box
        flag = 0
        if abs(p0.x() - self.xmin) < radius: flag = 1
        if abs(p0.x() - self.xmax) < radius: flag = 1
        if abs(p0.y() - self.ymin) < radius: flag = 1
        if abs(p0.y() - self.ymax) < radius: flag = 1
        
        # if near edge, manually check all objects, else just do the regular kd-tree search
        if flag and checkEdges:
            
            # do it manually by sorting and finding only objects with seps less than radius
            neighbor_list, seps_sorted = p0.sort(self.objects, [], [], maxDist=radius, projected=True)  
            s = set(neighbor_list)
        else:
            # do the recursive range search
            self._rangeSearch(self._root, rect, s, p0, radius)        
            
        # make a copy of the set
        scopy = copy.deepcopy(s)
                
        # run the function on the set, if specified
        if functionToRun is not None:
            for x in s:
                functionToRun(x)
       
        return scopy
        
    def _checkNearestSearch(self, n, p0, closest, depth):
        """
        @brief: check that point object closest is in fact the nearest neighbor of point p
        @param n: Node object that is head of tree
        @param p: point object that we are looking for the neighbor of
        @param closest: the believed to be closest point
        @param depth: current depth in the subtree
        """
        
        dist = p.distanceTo(n._element)    # dist from p to curr point
        min = p.distanceTo(closest)  # minimum distance
        
        newClosest = closest.copy()
        thisp = n._element # point corresponding to this node
        
        # check if curr dist is smaller than minimum
        if (dist < min):  
            newClosest = thisp
        
        if (n._lb is not None): 
            # only check left subtree if distance from p to the 
            # bounding box of left child is smaller than min
            if (n._lb._rect.distanceTo(p0) < min):
                newClosest = self._checkNearestSearch(n._lb, p0, newClosest, depth+1)
    
        if (n._rt is not None):
            # only check right subtree if distance from p to the 
            # bounding box of right child is smaller than min
            if (n._rt._rect.distanceTo(p0) < min):
                newClosest = self._checkNearestSearch(n._rt, p0, newClosest, depth+1)
        
        return newClosest

    def _nearestSearch(self, n, p0, closest, depth):
        """
        @brief: estimate the nearest neighbor of Point p by traversing the path taken
                when inserting p into the kdtree
        @param n: Node object that is head of tree
        @param p: point object that we are looking for the neighbor of
        @param closest: the believed to be closest point object
        @param depth: current depth in the subtree
        """
    
        # return if at the end of the tree
        if (n is None): 
            return closest
        
        newClosest = closest.copy()
        thisp = n._element # point corresponding to this node
        
        min = p0.distanceTo(closest) # minimum distance
        dist = p0.distanceTo(n._element) 
        
        # compare Nodes by x coord if at even depth 
        if ((depth % 2) == 0):
            
            # go to the left subtree
            if (p0.x() < thisp.x()):
                
                # check if distance is smaller and continue recursion
                if (dist < min): 
                    newClosest = self._nearestSearch(n._lb, p0, thisp, depth+1) 
                else:
                    newClosest = self._nearestSearch(n._lb, p0, closest, depth+1)
            
            # go to the right subtree
            else:
                
                # check if distance is smaller and continue recursion
                if (dist < min): 
                    newClosest = self._nearestSearch(n._rt, p0, thisp, depth+1) 
                else:
                    newClosest = self._nearestSearch(n._rt, p0, closest, depth+1)
        
        # compare Nodes by y coord
        else:
            
            # go to the left subtree
            if (p0.y() < thisp.y()):
            
                # check if distance is smaller and continue recursion
                if (dist < min): 
                    newClosest = self._nearestSearch(n._lb, p0, thisp, depth+1) 
                else:
                    newClosest = self._nearestSearch(n._lb, p0, closest, depth+1)  
            # go the right subtree
            else:
                
                # check if distance is smaller and continue recursion
                if (dist < min): 
                    newClosest = self._nearestSearch(n._lb, p0, thisp, depth+1)
                else:
                    newClosest = self._nearestSearch(n._lb, p0, closest, depth+1)
                    
        return newClosest
        
    def nearest(self, obj):
        """
        @brief: return a nearest neighbor in the set to obj; None if set is empty
        """
        
        # return None if empty
        if (self.isEmpty()):
            return None
        
        # make the point object with necessary functions
        p0 = self.pointClass(obj.__dict__.copy())
        
        # estimate nearest neighbor with method #1
        closest = self._nearestSearch(self._root, p0, self._root._element, 0)
        
        # confirm nearest neighbor
        closest = self._checkNearestSearch(self._root, p0, closest, 0)
        
        return closest

        
    
