'''
Created on December 2nd 2013

Objects that represent each TE annotation
Stores information regarding:
- start and end coordinates
- supporting information (not implemented yet)

@author: cjustin
'''

class TransEle:
    
    _transID = ""
    _chr = 0
    _start = 0
    _end = 0
    
    def __init__(self, transID, chr, start, end):
        self._transID = transID
        self._chr = chr
        self._start = start
        self._end = end
    
    def getTransID(self):
        return self._transID    
    
    def getStart(self):
        return self._start
    
    def getEnd(self):
        return self._end