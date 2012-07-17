__author__ = "Adam Traub"
__date__="1/21/2011"

"""
The Singly Linked List is a well established Data Structure in Computing.
It's not very efficient, but it's not contiguous in memory unlike arrays;
which means that the efficiency hit that arrays takewhen new memory needs
to be allocated can be avoided.  Adding an item to the front of the List
is an O(1) operation. Appending an item is O(n).
"""

class LinkedList:
    '''Linked List'''
    class __Node:
        '''private node object.  Node's consist of an element and a pointer to the next Node'''
        def __init__(self, element=None):
            self.element = element
            self.next = None

        def __str__(self):
            '''string representation of Node'''
            return str(self.element)

        def hasNext(self):
            '''returns true if the Node points to another Node'''
            return self.next != None

        def getNext(self):
            '''get next Node'''
            return self.next

        def getElement(self):
            '''get the element of the currentNode'''
            return self.element

        def setNext(self,nextItem):
            '''set the next Node of the currentNode'''
            self.next = nextItem

        def setElement(self, element):
            '''set the element of the current Node'''
            self.element = element

    def __init__(self):
        '''Creates the LinkedList'''
        self.first = LinkedList.__Node()
        self.length = 0

    def __len__(self):
        '''returns the length of the list O(1)'''
        return self.length

    def isEmpty(self):
        '''returns true if the list is empty'''
        return self.length == 0

    def append(self, element):
        '''Add element to last position of the list'''
        currentNode = self.first
        while currentNode.hasNext():
            currentNode = currentNode.getNext()
        currentNode.setNext(LinkedList.__Node(element))
        self.length+=1

    def popFront(self):
        '''Optimized method to pop front element from list O(1)'''
        if self.length > 0:
            val = self.first.getNext().getElement()
            self.first.setNext(self.first.getNext().getNext())
            self.length -= 1
            return val
        else:
            raise IndexError("Cannot pop from an empty List")

    def addToFront(self,element):
        '''Optimized method to add an element to the front of the list O(1)'''
        toAdd = LinkedList.__Node(element)
        toAdd.setNext(self.first.getNext())
        self.first.setNext(toAdd)
        self.length += 1

    def __getitem__(self, index):
        '''Allows for indexing, index must be an integer'''
        #accounts for slicing
        if type(index) == slice:
            return self.__sliceList(index)

        return self.__getNodeAtPosition(self.__checkIndex(index)).getElement()

    def __add__(self, other):
        retList = LinkedList()
        for item in self:
            retList.append(item)
        for item in other:
            retList.append(item)
        return retList
                    
    def pop(self, index =None):
        '''Removes and returns an element in the list. last element by default.'''
        if self.length == 0:
            raise IndexError("pop from empty list")
        
        #utilize default parameter    
        if index ==None: 
            index = self.length-1

        previous = self.__getNodeAtPosition(self.__checkIndex(index)-1)
        toRemove = previous.getNext()
        afterNext = None
        if toRemove.hasNext():
            afterNext = toRemove.getNext()
        previous.setNext(afterNext)
        self.length-=1
        return toRemove.getElement()

    def where(self, element):
        """ return the index of where the element exists in the list"""
        index = 0
        for item in self:
            if item == element:
                return index
            index += 1
            
        return None
        
    def __setitem__(self, index, element):
        '''Sets the item at a given index to a new element.'''
        self.__getNodeAtPosition(self.__checkIndex(index)).setElement(element)

    def __str__(self):
        '''returns a string representation of the list'''
        if self.length == 0:    
            return '[]'
        retString = "["
        currentElement = self.first.getNext()
        for i in range(self.length):
            retString += str(currentElement) +", "
            currentElement = currentElement.getNext()
        
        return retString[:-2] + ']'

    def insert(self, element, index):
        '''inserts an element to the given index'''
        previous = self.__getNodeAtPosition(self.__checkIndex(index)-1)
        toMove = previous.getNext()
        previous.setNext(LinkedList.__Node(element))
        previous=previous.getNext()
        previous.setNext(toMove)
        self.length+=1


    #Private functions
    def __sliceList(self,theSlice):
        '''(Private) function to handle slicing.  Returns a Linked List'''
        retList = LinkedList()
        
        #Following conditions handles the x[start:stop:step] notation
        step  = self.__determineStartStopStep(theSlice.step,1,1,1)
        start = self.__determineStartStopStep(theSlice.start,step,0,self.length-1)
        stop  = self.__determineStartStopStep(theSlice.stop,step,self.length,-1)

        for eachItem in range(start,stop,step):
            retList.append(self.__getitem__(eachItem))
        return retList

    def __determineStartStopStep(self, valToCheck, step,
                             positiveStep, negativeStep):
        '''private function to reduce repeated code in determining slicing handling'''
        if valToCheck == None:
            if step > 0:
                return positiveStep
            else:
                return negativeStep
        else:
            return valToCheck
        
    def __getNodeAtPosition(self, index):
        '''(Private) Gets a Node at a given index'''
        currentNode = self.first
        for i in range(index+1):#Adds 1 to account for initial Node (sentinel)
            currentNode = currentNode.getNext()
        return currentNode

    def __checkIndex(self, index):
        '''(Private) check if the index is an acceptable value.  Index only changes if negative'''
        if type(index) != int:
            raise TypeError("Index must be an integer or a slice not a "+str(type(index)).split("'")[1])

        #handles negative indices.
        if index < 0:
            index += self.length

        #If the index is out of bounds
        if index >= self.length or index < 0:
            raise IndexError("Index out of bounds")
    
        return index
