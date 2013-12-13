'''
Created on December 12th 2013

Objects that represent each exon annotation
Stores information regarding:
- exon coordinates
- coverage counts for the alleles

@author: cjustin
'''

class ExonInfo:

    #start and end locations are used to define exon
    _exonStart = 0
    _exonEnd = 0
    
    #m1 = maternal strain1
    _exonCov_m1 = 0
    _exonCov_p2 = 0
    _exonCov_p1 = 0
    _exonCov_m2 = 0

    """
    Constructor
    """
    def __init__(self, start, end):
        self._exonStart = int(start)
        self._exonEnd = int(end)
        self._exonCov_m1 = 0.0
        self._exonCov_p2 = 0.0
        self._exonCov_p1 = 0.0
        self._exonCov_m2 = 0.0
    
    def getStart(self):
        return self._exonStart

    def getEnd(self):
        return self._exonEnd   
    """
    Returns coverage of gene for exon for m1
    """
    def setCoverage_m1(self, cov):
         self._exonCov_m1 = float(cov)
    
    """
    Returns coverage of gene for exon for p2
    """
    def setCoverage_p2(self, cov):
        self._exonCov_p2 = float(cov)

    """
    Returns coverage of gene for exon for p1
    """
    def setCoverage_p1(self, cov):
        self._exonCov_p1 = float(cov)

    """
    Returns coverage of gene for exon for m2
    """
    def setCoverage_m2(self, cov):
        self._exonCov_m2 = float(cov)
    
    
    """
    Returns coverage of gene for exon for m1
    """
    def getCoverage_m1(self):
        return self._exonCov_m1
    
    """
    Returns coverage of gene for exon for p2
    """
    def getCoverage_p2(self):
        return self._exonCov_p2

    """
    Returns coverage of gene for exon for p1
    """
    def getCoverage_p1(self):
        return self._exonCov_p1

    """
    Returns coverage of gene for exon for m2
    """
    def getCoverage_m2(self):
        return self._exonCov_m2
    
