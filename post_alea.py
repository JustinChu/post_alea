'''
Created on December 2nd 2013

inputs:
-g gene_annotations_bed
-1 strain1_TE_annotations_bed
-2 strain2_TE_annotations_bed

#Originally from same set of reads
--m1 strain1_maternal_bedGraph file
--p2 strain2_paternal_bedGraph file

#Originally from same set of reads
--p1 strain1_paternal_bedGraph file
--m2 strain2_maternal_bedGraph file

@author: cjustin
'''

from gene_info import GeneInfo
from optparse import OptionParser
import os
import re
# import math
import track

class PostAlea:
    
    # header values in array for each of changing order
    # coverages m1:p2:p1:m2
    _headerFields = ['event_ID', 'TE_name', 'UCSC_gene', 'TE_location', 'gene_location',
                     'distance', 'ratio', 'coverage_total', 'coverages', 'exon_info']
    
    # fileNameStorage
    _te1FileName = ""
    _te2FileName = ""
    
    _geneList = {}
    _geneNames = {}
    
    # filtering and search variables
    _coverageThreshold = 0
    _windowSize = 0
    
    """
    Constructor
    """
    def __init__(self, gene, te1, te2, windowSize, coverageThreshold):
        self._te1FileName = te1
        self._te2FileName = te2
        self._windowSize = windowSize
        self._coverageThreshold = coverageThreshold
        self._geneNames = self._getGeneNames(gene)
    
    """
    Extracts the USCS gene IDs for each refseqID for later lookup
    Uses a refFlat file of refseq genes
    """
    def _getGeneNames(self, refFlat):
        
        geneNames = {}
        
        fh = open(refFlat)
        
        headerArray = fh.readline().split("\t")
        
        # get header indexes
        indexRefSeqName = -1
        indexGeneName2 = -1
        
        for i in range(len(headerArray)):
            if headerArray[i] == "name":
                indexRefSeqName = i
            if headerArray[i] == "name2":
                indexGeneName2 = i
        
        for i in fh:
            lineHeader = i.split("\t")
            geneNames[lineHeader[indexRefSeqName]] = lineHeader[indexGeneName2]
        
        return geneNames
    
    """
    Starts the analysis
    """                    
    def run(self, m1, p2, p1, m2, outputPrefix):
        # populate list of genes that show any uniqueness with regards to allelic ratios
        
        # make filehandles for each file
        fh_m1 = open(m1)

        # create GeneInfo objects
        for i in fh_m1:            
            # break line into array
            lineArray = i.strip().split("\t")
            
            currentGeneID = lineArray[3]
            currentChr = lineArray[0]

            # add gene to list if it does not yet exist            
            if not self._geneList.has_key(currentGeneID):
                geneInfoObj = GeneInfo(currentGeneID, currentChr)
                self._geneList[currentGeneID] = geneInfoObj

            # populate exon list
            currentStart = lineArray[1]
            currentEnd = lineArray[2] 
            self._geneList[currentGeneID].addUpdateExon(currentStart, currentEnd,
                                                        lineArray[5], 0, 0, 0)
        # close file handle
        fh_m1.close()
        
        fh_p2 = open(p2)
        for i in fh_p2:
            # break line into array
            lineArray = i.strip().split("	")
            
            currentGeneID = lineArray[3]
            currentChr = lineArray[0]

            # add gene to list if it does not yet exist            
            if not self._geneList.has_key(currentGeneID):
                geneInfoObj = GeneInfo(currentGeneID, currentChr)
                self._geneList[currentGeneID] = geneInfoObj

            # populate exon list
            currentStart = lineArray[1]
            currentEnd = lineArray[2] 
            self._geneList[currentGeneID].addUpdateExon(currentStart, currentEnd,
                                                        0, lineArray[5], 0, 0)
        fh_p2.close()
 
        fh_p1 = open(p2)
        for i in fh_p1:
            # break line into array
            lineArray = i.strip().split("\t")
            
            currentGeneID = lineArray[3]
            currentChr = lineArray[0]

            # add gene to list if it does not yet exist            
            if not self._geneList.has_key(currentGeneID):
                geneInfoObj = GeneInfo(currentGeneID, currentChr)
                self._geneList[currentGeneID] = geneInfoObj

            # populate exon list
            currentStart = lineArray[1]
            currentEnd = lineArray[2] 
            self._geneList[currentGeneID].addUpdateExon(currentStart, currentEnd,
                                                        0, 0, lineArray[5], 0)
        fh_p1.close()
        
        fh_m2 = open(m2)
        for i in fh_m2:
            # break line into array
            lineArray = i.strip().split("\t")
            
            currentGeneID = lineArray[3]
            currentChr = lineArray[0]

            # add gene to list if it does not yet exist            
            if not self._geneList.has_key(currentGeneID):
                geneInfoObj = GeneInfo(currentGeneID, currentChr)
                self._geneList[currentGeneID] = geneInfoObj

            # populate exon list
            currentStart = lineArray[1]
            currentEnd = lineArray[2] 
            self._geneList[currentGeneID].addUpdateExon(currentStart, currentEnd,
                                                        0, 0, 0, lineArray[5])
        fh_m2.close()

        strain1_TEs = track.load(self._te1FileName)
        strain2_TEs = track.load(self._te2FileName)
        
        # stores all events to finally output
        # eventID is the geneName_TEtype_TEposition
        # eventID->fields->values
        events = {}
        
        # eventsFiltered = {}

        # genes with any evidence of allelic skew
        # look for TEs
        for i in self._geneList.keys():
            gene = self._geneList[i]
            coverageTotal = gene.getCoverage()
            if coverageTotal > self._coverageThreshold:
                ratio = gene.getTotalAllelicRatio()
                coverages = gene.getSummarizedCoveragesStr()
                
                # for each gene look with window size for TEs
                candidateTEs1 = strain1_TEs.read({'chr':gene.getChr(),
                                                  'start':(gene.getStart() - self._windowSize),
                                                  'end':(gene.getEnd() + self._windowSize)})
                
                # 'event_ID', 'TE_location', 'gene_location', 'distance', 'ratio', 'coverage_total', 'coverages'
                for event in candidateTEs1:
                    eventHash = {}
                    eventStart = event[0]
                    eventEnd = event[1]
                    eventType = event[2]
                    eventID = gene.getGeneID() + "_" + str(eventStart) + "_" + eventType
                    
                    eventHash["event_ID"] = eventID
                    eventHash["TE_name"] = eventType
                    eventHash["TE_location"] = gene.getChr() + ":" + str(eventStart) + "-" + str(eventEnd)
                    eventHash["gene_location"] = gene.getChr() + ":" + str(gene.getStart()) + "-" + str(gene.getEnd())
                    eventHash["distance"] = eventStart - gene.getStart()
                    eventHash["ratio"] = ratio
                    eventHash["coverage_total"] = coverageTotal
                    eventHash["coverages"] = coverages
                    eventHash["UCSC_gene"] = self._geneNames[gene.getGeneID()]
                    eventHash["exon_info"] = self._geneNames[gene.getExonInfo()]
                    
                    events[eventID] = eventHash
                
                candidateTEs2 = strain2_TEs.read({'chr':gene.getChr(),
                                                  'start':(gene.getStart() - self._windowSize),
                                                  'end':(gene.getEnd() + self._windowSize)})
                for event in candidateTEs2:
                    eventHash = {}
                    eventStart = event[0]
                    eventEnd = event[1]
                    eventType = event[2]
                    eventID = gene.getGeneID() + "_" + str(eventStart) + "_" + eventType
                    
                    eventHash["event_ID"] = eventID
                    eventHash["TE_name"] = eventType
                    eventHash["TE_location"] = gene.getChr() + ":" + str(eventStart) + "-" + str(eventEnd)
                    eventHash["gene_location"] = gene.getChr() + ":" + str(gene.getStart()) + "-" + str(gene.getEnd())
                    eventHash["distance"] = gene.getStart() - eventStart
                    eventHash["ratio"] = ratio
                    eventHash["coverage_total"] = coverageTotal
                    eventHash["coverages"] = coverages
                    eventHash["UCSC_gene"] = self._geneNames[gene.getGeneID()]
                    eventHash["exon_info"] = self._geneNames[gene.getExonInfo()]
                    
                    events[eventID] = eventHash                
        
        # output file
        out_fh = open(outputPrefix + '.tsv', 'w')
        out_fh.write(self._getHeader())
        
        for event in events.keys():
            out_fh.write(self._getLine(events[event]))
               
        # close file handles
        out_fh.close()
        
        # output file for filtered events
        # out_filtered_fh = open(outputPrefix + 'filtered.tsv', 'w')
        # out_filtered_fh.close()
    
    """
    Helper method
    Returns header as a string for printing
    """            
    def _getHeader(self):
        headerStr = ""
        for field in self._headerFields:
            headerStr += field + "\t"
        return headerStr.strip() + "\n"
    
    """
    Helper method
    params eventHash - Dict a single event. Each key should 
                       be one of the header array values.
    Returns eventHash contents as a tab delimited string
    """
    def _getLine(self, eventHash):
        outStr = ""
        for field in self._headerFields:
            outStr += str(eventHash[field]) + "\t"
        return outStr.strip() + "\n"

if __name__ == '__main__':
    parser = OptionParser()
    parser.add_option("-g", "--gene", dest="gene",
                      help="refFlat file of gene annotations (Required)", metavar="GENEBED")
    parser.add_option("-1", "--strain1_TE", dest="te1", metavar="TE1",
                      help="bed file of Transposible Elements (TE) for strain1 (Required)")
    parser.add_option("-2", "--strain2_TE", dest="te2", metavar="TE2",
                      help="bed file of Transposible Elements (TE) for strain2 (Required)")
    parser.add_option("--m1", dest="m1", metavar="M1",
                      help="coverage file for strain1 where it was the maternal strain (Required)")
    parser.add_option("--p2", dest="p2", metavar="P2",
                      help="coverage file for strain2 where it was the paternal strain (Required)")
    parser.add_option("--p1", dest="p1", metavar="P1",
                      help="coverage file for strain1 where it was the paternal strain (Required)")
    parser.add_option("--m2", dest="m2", metavar="M2",
                      help="coverage file for strain2 where it was the maternal strain (Required)")
    parser.add_option("-c", "--coverage_threhold", dest="coverageThreshold",
                      metavar="COVERAGE_THRESHOLD", default=50,
                      help="coverage threshold needed evaluate event")
    parser.add_option("-w", "--window_size", dest="windowSize",
                      metavar="WINDOW_SIZE", default=20000,
                      help="distance from gene to consider transposable elements")
    parser.add_option("-o", "--output_prefix", dest="outputPrefix",
                      metavar="WINDOW_SIZE", default="output",
                      help="output Prefix to output files in")
    
    (options, args) = parser.parse_args()
    
    if (options.gene and options.te1 and options.te2 
        and options.m1 and options.p2 and options.p1 and options.m2):
        runner = PostAlea(options.gene, options.te1, options.te2,
                          options.windowSize, options.coverageThreshold)
        runner.run(options.m1, options.p2, options.p1, options.m2, options.outputPrefix)
    else:
        print 'ERROR: Missing Required Options. Use -h for help'

