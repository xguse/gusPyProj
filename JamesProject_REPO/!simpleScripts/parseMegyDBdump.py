'''
This script is intended to parse the following file format, as supplied by Karyn Megy from VectorBase:

******
SLICE: -supercont3.1-
GENE: CPIJ020333 - K.supercont3.1: pos.15291-15363 (-1)
	TRANSCRIPT: CPIJ020333-RA - K.supercont3.1: pos.15291-15363 (-1)
		EXON: E073991 - K.supercont3.1: pos.15291-15363 (-1)
******
SLICE: -supercont3.1-
GENE: CPIJ000001 - K.supercont3.1: pos.20593-42045 (-1)
	TRANSCRIPT: CPIJ000001-RA - K.supercont3.1: pos.20593-42045 (-1)
		EXON: E014733 - K.supercont3.1: pos.41821-42045 (-1)
		  EXON CODING REGION: 41821-42045
		EXON: E014734 - K.supercont3.1: pos.41621-41755 (-1)
		  EXON CODING REGION: 41621-41755
		EXON: E014735 - K.supercont3.1: pos.41382-41445 (-1)
		  EXON CODING REGION: 41382-41445
		EXON: E014736 - K.supercont3.1: pos.20593-20735 (-1)
		  EXON CODING REGION: 20593-20735
******
SLICE: -supercont3.1-
GENE: CPIJ020335 - K.supercont3.1: pos.77166-77239 (1)
	TRANSCRIPT: CPIJ020335-RA - K.supercont3.1: pos.77166-77239 (1)
		EXON: E073992 - K.supercont3.1: pos.77166-77239 (1)
******
SLICE: -supercont3.1-
GENE: CPIJ020336 - K.supercont3.1: pos.93017-93084 (-1)
	TRANSCRIPT: CPIJ020336-RA - K.supercont3.1: pos.93017-93084 (-1)
		EXON: E073993 - K.supercont3.1: pos.93017-93084 (-1)
'''
import string

class megyFeatures:
    """
    Takes   -> a single gene feature str from above file type.
    
    Does    -> parses string into features
            -> provides methods to output a genes exon coords as my coding bounds scripts expect
    """
    # Attribs
    formatedRec = None
    _workingRec = None
    _parsedRec  = {}
    
    def __init__(self, geneFeatureStr):
        self._workingRec = geneFeatureStr
        self._parseFtrStr()
        pass
    
    def _parseFtrStr(self):
        # split on '\n'
        workingRec = self._workingRec.split('\n')
        parsedRec  = {'SLICE':'',
                      'GENE':[],
                      'TRANSCRIPTS':[],
                      'EXONS':[]
                      }
        i = 0
        numPfLines = len(workingRec)
        while i < numPfLines:
            if workingRec[i].startswith('SLICE'):
                parsedRec['SLICE'] = workingRec[i].split('-')[1]
                i+=1
            elif workingRec[i].startswith('GENE'):
                geneID     = workingRec[i][6:16]
                gfields     = workingRec[i].split(' ')
                start,stop = (gfields[4].replace('.','-').split('-')[1:3])   #parses 'pos.20593-42045' to (20593,42045) which is unpacked into respective vars
                strand     = gfields[5].replace(')','(').split('(')[1]
                parsedRec['GENE'].extend([geneID,start,stop,strand])
                i+=1
            elif workingRec[i].strip().startswith('TRANSCRIPT'):
                transcriptID = workingRec[i][14:26]
                tfields       = workingRec[i].split(' ')
                start,stop   = (tfields[4].replace('.','-').split('-')[1:3])   #parses 'pos.20593-42045' to (20593,42045) which is unpacked into respective vars
                parsedRec['TRANSCRIPTS'].extend([transcriptID,start,stop])
                i+=1
            elif workingRec[i].strip().startswith('EXON'):
                efields   = workingRec[i].split(' ')
                exonID       = efields[1]
                start,stop   = (efields[4].replace('.','-').split('-')[1:3])   #parses 'pos.20593-42045' to (20593,42045) which is unpacked into respective vars
                cStart,cStop = '',''
                if workingRec[i+1].strip().startswith('EXON CODING'):
                    i+=1
                    cStart,cStop = (workingRec[i].split(' ')[5].split('-'))
                parsedRec['EXONS'].append([exonID,start,stop,cStart,cStop])
                i+=1
            else:
                i+=1
        
        # sort Exons so they are all ordered such that
        # the exon starting with the smallest coord is
        # listed first
        parsedRec['EXONS'].sort(key=lambda x: x[1])
        
        # format a multi-line string to represent the final record
        finalRec = []
        for exon in parsedRec['EXONS']:
            finalRec.append('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (parsedRec['GENE'][0],parsedRec['SLICE'],'--',parsedRec['GENE'][3],parsedRec['TRANSCRIPTS'][0],exon[3],exon[4],exon[1],exon[2],exon[0]))
        
        # assign rec to class var
        self.formatedRec = ''.join(finalRec)
        ##print self.formatedRec
    

    
iFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/Culex_CpipJ1-2_GeneTranscrExon.txt'
oFile = '/Users/biggus/Documents/James/Data/Tu_miRNA/Cq_CpipJ1-2_GeneTranscrExon_112408.txt'
sep   = '******'

# read in file data, format as one long string, then chop it up into a list of 'gene' records
iFile = open(iFile,'rU').readlines()
iFile = ''.join(iFile).split(sep)


outPut = ['#geneID\tchrom\tbiotype\tstrand\ttransID\tcodingStart\tcodingEnd\texonStart\texonStop\texonID\n']

# main loop
c=1
for gene in iFile:
    geneFtrs = megyFeatures(gene)
    formatedFtrs = geneFtrs.formatedRec
    outPut.append(formatedFtrs)
    print c
    c+=1
    

# write out
oFile = open(oFile,'w')
oFile.writelines(outPut)
oFile.close()




