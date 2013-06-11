"""
 random_sampling.py
 implements random sampling with and without replacement, as described 
 here: http://stackoverflow.com/questions/2140787/ \
 select-random-k-elements-from-a-list-whose-elements-have-weights/2149533#2149533
 
 author: Nick Hand
 contact: nhand@berkeley.edu
 creation date: 06/11/2013
"""
import random

class Node:
    """
    A class to implement a node in the heap which as a weight, value and 
    total weight. The total weight, self.tw, is self.w plus the weight of 
    any children
    """
    __slots__ = ['w', 'v', 'tw']

    def __init__(self, w, v, tw):
        self.w, self.v, self.tw = w, v, tw
#endclass Node

#-------------------------------------------------------------------------------
def rws_heap(items):
    """
    Initialize a binary heap with dict items. The heap is a binary heap that lives 
    in an array. It has a Node for each pair in 'items.'  
    
    Parameters
    ----------
    items : list of tuples
        The items to store. The tuples should contain the (value, weight)
        of each item you wish to sample from
        
    Notes
    -----
    h[1] is the root. Each other Node h[i] has a parent at h[i>>1]. 
    Each node has up to 2 children, h[i<<1] and h[(i<<1)+1]. To get this nice 
    simple arithmetic, we have to leave h[0] vacant.
    """
    h = [None]                          # leave h[0] vacant
    for v, w in items:
        h.append(Node(w, v, w))

    for i in range(len(h) - 1, 1, -1):  # total up the tws
        h[i>>1].tw += h[i].tw           # add h[i]'s total to its parent

    return h
#end rws_heap

#-------------------------------------------------------------------------------
def rws_heap_pop(h):
    """
    Pop an item from the heap
    """
    gas = h[1].tw * random.random()      # start with a random amount of gas

    i = 1                  # start driving at the root
    while gas > h[i].w:    # while we have enough gas to get past node i:  
        gas -= h[i].w      #      drive past node i
        i <<= 1            #      move to first child
        if gas > h[i].tw:  #      if we have enough gas:
            gas -= h[i].tw #          drive past first child and descendants
            i += 1         #          move to second child

    w = h[i].w             # out of gas! h[i] is the selected node
    v = h[i].v

    h[i].w = 0             # make sure this node isn't chosen again
    while i:               # fix up total weights
        h[i].tw -= w
        i >>= 1

    return v
#end rws_heap_pop

#-------------------------------------------------------------------------------
def random_weighted_sample_no_replacement(items, N):
    """
    Randomly choose N items by the weights associated, with no replacement 
        
    Parameters
    ----------
    items : list of tuples
        The items to store. The tuples should contain the (value, weight)
        of each item you wish to sample from 
    N : int
        the number of objects to choose
    """
    heap = rws_heap(items)               # just make a heap...
    objs = []
    for i in xrange(N):
        objs.append(rws_heap_pop(heap))  # and pop n items off it

    return objs
#end random_weighted_sample_no_replacement

#-------------------------------------------------------------------------------
def random_weighted_sample_with_replacement(items, N):
    """
    Randomly choose N items by the weights associated, with replacement
    
    Parameters
    ----------
    items : list of tuples
        The items to store. Tuples should contain the (value, weight)
        of each item you wish to sample from 
    N : int
        the number of objects to choose
    """
    total = float(sum(w for v, w in items))
    i = 0
    v, w = items[0]
    objs = []
    while n:
        
        x = total * (1 - random.random() ** (1.0 / n))
        total -= x
        while x > w:
            x -= w
            i += 1
            v, w = items[i]
        w -= x
        objs.append(v)
        n -= 1

    return objs
#end random_weighted_sample_with_replacement

#-------------------------------------------------------------------------------    

