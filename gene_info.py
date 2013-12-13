'''
Created on December 2nd 2013

Objects that represent each gene annotation
Stores information regarding:
- start and end coordinates
- exon coordinates
- coverage counts for both alleles

@author: cjustin
'''

from exon_info import ExonInfo

class GeneInfo:
    
    _geneID = ""
    _chr = 0
    _totalCov = 0
    
    _start = 'Inf'
    _end = 0
    
    # Exon information
    _exonHash = {}
    
    """
    Constructor
    """
    def __init__(self, geneID, chr):
        self._geneID = geneID
        self._chr = chr
        self._totalBases = 0
        self._totalCov = 0
        self._exonHash = {}
        self._start = 'Inf'
        self._end = 0

    """
    Adds or updates a single exon geneInfo object
    Exon information consists of coverage and location information
    
    Coverages in each update will be overwritten if the new coverage is
    greater than 0
    """         
    def addUpdateExon(self, start, end, m1, p2, p1, m2):
        
        exonID = start + "_" + end
        
        if not self._exonHash.has_key(exonID):
            self._exonHash[exonID] = ExonInfo(start, end)
            if int(start) < self._start:
                self._start = int(start)
            if int(end) > self._end:
                self._end = int(end)
            self._totalBases += int(end) - int(start)
        
        if m1 != 0 :
            self._exonHash[exonID].setCoverage_m1(m1)
        if p2 != 0 :
            self._exonHash[exonID].setCoverage_p2(p2)
        if p1 != 0 :
            self._exonHash[exonID].setCoverage_p1(p1)
        if m2 != 0 :
            self._exonHash[exonID].setCoverage_m2(m2)
                
    def getGeneID(self):
        return self._geneID
    
    def getChr(self):
        return self._chr
    
    def getStart(self):
        return self._start
    
    def getEnd(self):
        return self._end
    
    
    """
    Returns normalized coverage of gene for all exons
    """
    def getCoverage(self):
        cov = 0
        
        for i in self._exonHash.keys():
            # proportion of gene exon is responsible for
            geneProportion = (self._exonHash[i].getEnd() - self._exonHash[i].getStart()) / self._totalBases
            cov += (self._exonHash[i].getCoverage_m1()) * geneProportion
        return cov  
    
    """
    Returns normalized coverage of gene for all exons for m1
    """
    def getCoverage_m1(self):
        cov = 0
        
        for i in self._exonHash.keys():
            # proportion of gene exon is responsible for
            geneProportion = (self._exonHash[i].getEnd() - self._exonHash[i].getStart()) / self._totalBases
            cov += self._exonHash[i].getCoverage_m1() * geneProportion
        return cov
    
    """
    Returns normalized coverage of gene for all exons for p2
    """
    def getCoverage_p2(self):
        cov = 0
        
        for i in self._exonHash.keys():
            # proportion of gene exon is responsible for
            geneProportion = (self._exonHash[i].getEnd() - self._exonHash[i].getStart()) / self._totalBases
            cov += self._exonHash[i].getCoverage_p2() * geneProportion
        return cov

    """
    Returns normalized coverage of gene for all exons for p1
    """
    def getCoverage_p1(self):
        cov = 0
        
        for i in self._exonHash.keys():
            # proportion of gene exon is responsible for
            geneProportion = (self._exonHash[i].getEnd() - 
                              self._exonHash[i].getStart()) / self._totalBases
            cov += self._exonHash[i].getCoverage_p1() * geneProportion
        return cov

    """
    Returns normalized coverage of gene for all exons for m2
    """
    def getCoverage_m2(self):
        cov = 0
        
        for i in self._exonHash.keys():
            # proportion of gene exon is responsible for
            geneProportion = (self._exonHash[i].getEnd() - 
                              self._exonHash[i].getStart()) / self._totalBases
            cov += self._exonHash[i].getCoverage_m2() * geneProportion
        return cov
    
    def getSummarizedCoveragesStr(self):
        return str(self.getCoverage_m1()) + ";" + str(self.getCoverage_p2()) \
            + ";" + str(self.getCoverage_p1()) + ";" + str(self.getCoverage_m2())
    
    """
    Returns normalized coverage of gene for all exons for strain 1
    """
    def getCoverageStrain1(self):
        cov = 0
        for i in self._exonHash.keys():
            # proportion of gene exon is responsible for
            geneProportion = (self._exonHash[i].getEnd() - 
                              self._exonHash[i].getStart()) / self._totalBases
            cov += (self._exonHash[i].getCoverage_m1() + 
                    self._exonHash[i].getCoverage_p1()) * geneProportion
        return cov    
    
    """
    Returns normalized coverage of gene for all exons for strain 2
    """
    def getCoverageStrain2(self):
        cov = 0
        for i in self._exonHash.keys():
            # proportion of gene exon is responsible for
            geneProportion = (self._exonHash[i].getEnd() - 
                              self._exonHash[i].getStart()) / self._totalBases
            cov += (self._exonHash[i].getCoverage_m2() + 
                    self._exonHash[i].getCoverage_p2()) * geneProportion
        return cov 
        
    """
    Computes Allelic Ratio
    
    Returns:
    >0.5 is biased to strain1
    <0.5 is biased to strain2    
    if -1 is returned then the gene no coverage to at all to compute ratio
    """
    def getTotalAllelicRatio(self):
        covStrain1 = self.getCoverageStrain1()
        covStrain2 = self.getCoverageStrain2()
        
        if covStrain1 == 0 and covStrain2 == 0:
            return -1
        else:
            return covStrain1 / (covStrain1 + covStrain2)
        
    """
    Computes Allelic Ratio of maternal/paternal allele
    
    Returns:
    >0.5 is biased to maternal strains
    <0.5 is biased to paternal strains    
    if -1 is returned then the gene no coverage to at all to compute ratio
    """
    def getMatPatRatio(self):
        covStrain1 = 0
        covStrain2 = 0
        
        for i in self._exonHash.keys():
            # proportion of gene exon is responsible for
            geneProportion = (self._exonHash[i].getEnd() - 
                              self._exonHash[i].getStart()) / self._totalBases
            covStrain1 += (self._exonHash[i].getCoverage_m1() + 
                           self._exonHash[i].getCoverage_m2()) * geneProportion
            covStrain2 += (self._exonHash[i].getCoverage_p1() + 
                           self._exonHash[i].getCoverage_p2()) * geneProportion
        
        if covStrain1 != 0 and covStrain2 != 0:
            return covStrain1 / (covStrain1 + covStrain2)
        else:
            return -1
        
        
        
