"""Motif database parser.
    
Contains parser for the Transfac database (supports up to version 11.3)

Implemented by Paul Rigor, Kenny Daily 2009, 2010
Institute of Genomics and Bioinformatics
University of California, Irvine

All software and related material contained 
herein can be downloaded freely for academic,
non-commercial, research use only.

For any other use, please contact:
Pierre Baldi at pfbaldi _at_ ics uci edu. 

Copyright 2009-2010

Contact: prigor _at_ ics uci edu, kdaily _at_ ics uci edu

"""

import os
import sys
from exceptions import NotImplementedError

############
### Base classes for Parsers, Controllers and Data containers
############

_transfac_fields = ['AC', 'AS', 'BA', 'BF', 'BS', 'CC', 'CO', 'DE', 'DT', 'ID', 
                    'NA', 'P0', 'RA', 'RL', 'RN', 'RT', 'RX', 'VV']

class BaseMotifParser(object):
    """ Base class for motif database parser.
    
    Important member variables:
    self.motifs will contain the dictionary of motifs with key-value pairs
    for their respective motif data
    """
    
    def __init__(self,filename,DEBUG=False):
        self.motifs = BaseDataDict({})
        self.filename = filename
        self.rawinputlines = []
        self.DEBUG = DEBUG
        
    def getData(self):
        """Returns the motif dictionary.
        
        """
        
        return self.motifs
    
    def loadFile(self,filename=None):
        """Load the input file, if present.

        """
        
        if not filename:
            filename = self.filename
        try:
            ifile = open(filename,"r")
            self.rawinputlines = [line.strip() for line in ifile.readlines()]
        except Exception,e:
            raise Exception("Unable to read '%s': %s" % (filename,e.message))
        return self.rawinputlines
    
    def generateData(self):
        """ Need to overload this in the sub-class.
        
        """
        
        raise NotImplementedError



class BaseMotifController(object):
    def parse(self):
        raise NotImplementedError
    
    def getData(self):
        raise NotImplementedError

class BaseDataDict(dict):
    pass

class BaseMotifData(object):
    """Base class for accessing motif information after parsing step."""
    
    def __init__(self,data,DEBUG=False):
        self.data = BaseDataDict(data)
        self.DEBUG = DEBUG
    
    def getSectionData(self,motifkey,tag=None):
        """Obtain the information stored in 'tag' (return as list) using motifkey (accession id).
        """
        
        self.is_valid_motif(motifkey)
        data = []
        
        try:
            if tag!=None:
                data = self.data[motifkey][tag]
            else:
                data = self.data[motifkey]
        except:
            pass
        
        return data
    
    def setSectionData(self,motifkey,tag,data):
        """Sets the section data of a motif for a given section tag. 
        
        Be careful with this method.
        
        """
        
        self.is_valid_motif(motifkey)
        self.data[motifkey][tag]=data
        
    
    def getSectionList(self,motifkey):
        """Obtain the information stored in 'tag' (return as list) using motifkey (accession id).

        """
        self.is_valid_motif(motifkey)
        data = self.data[motifkey].keys()
        return data
        
    def is_valid_motif(self,motifkey):
        """Helper method to make ensure a valid motif accession.

        """
        
        if motifkey not in self.data.keys():
            raise KeyError("Cannot locate motif using accession id '%s'." % motifkey)
        return True            

    def is_valid_pwm(self,motifkey,pwm=None):
        """Checks whether the position weight matrix is a valid count, prob, freq matrix.

        """
        
        if pwm==None:
            pwm = self.getSectionData(motifkey,"pwm") # Fix this!
        if not pwm:
            validPWM=False
        else:
            hasDecimal="." in str(pwm)
            validPWM = False
            validProb = True
            validFreq = True
            if not hasDecimal:
                validPWM = True
            else:
                newPWM = []
                sums = []
                for row in pwm:
                    newRow = [float(col) for col in row]
                    theSum = sum(newRow)
                    if theSum != 1.0:
                        validProb = False  
                    newPWM.append(newRow)
                    sums.append(theSum)
                theSum = sums[0]
                for a_sum in sums[1:]:
                    if theSum!=a_sum:
                        validFreq = False
                if validProb or validFreq:
                    self.setSectionData(motifkey,tag,data)
                    validPWM=True
            #print "Freq: %s, Prob: %s" % (validFreq,validProb)
        return validPWM,validFreq,validProb                

    def uniquefy(self,alist):
        """Helper method to uniquefy lists.

        """
        
        u = {}
        for i in alist:
            u[i]=0
        
        return u.keys()            

    def __len__(self):
        """Returns the number of motifs stored.

        """
        
        return len(self.data.keys())

    def __getitem__(self,motifkey):
        """Obtain the annotations for a motif, usually a dictionary.

        """
        
        self.is_valid_motif(motifkey)
        return self.data[motifkey]
        
    def getShuffledPWM(self,motifkey,numlist):
        """Obtains the clean position weight matrix (without concensus column) and shuffles the rows
           to produce 'numlist' of shuffled matrices.
            
            returns a list of dictionaries:
                [ {'motifkey':
                    {'pwm':[
                            [], [], []
                            ]
                        } 
                    }
                ]    
        """
        pwmKey = "pwm"
        shuffledList = []
        pwm = self.getMotifPWMClean(motifkey)
        pwmLen = len(pwm)
        pwmIdxList = range(pwmLen)
        for counter in range(numlist):
            random.shuffle(pwmIdxList)
            # new shuffled pwm
            newPWM = []
            for idx in pwmIdxList:
                newPWM.append(pwm[idx])
            # append to the list of dictionaries
            shuffledList.append({motifkey:{pwmKey:newPWM}})
        return shuffledList
        

    def exportDataAsTransfac(self,motifname=None,species=None,pwm=None,stripSections=False):
        """ Export data to transfac format.
                motifname
                species
                pwm -- use this matrix
                
            Exports only AC,NA and P0 sections unless stripSections is set to False.
            
            Returns dataString
        """
        motifList = []
        if motifname==None:
            motifList = self.getMotifList(species)
        else:
            motifList.append(motifname)
        dataString = ""
        for motif in motifList:
            dataString += "//\n"
            dataString += "XX\n"
            dataString += "AC  %s\n" % motif
            dataString += "XX\n"
            dataString += "NA  %s\n" % self.getMotifGeneName(motif)
            dataString += "XX\n"
            if not stripSections:  #include all other sections as well
                sectionData = self.getSectionDataAsTransfac(motif)
                for tag in sectionData.keys():
                    if tag not in ['NA', 'pwm', 'AC']:
                        if hasattr(sectionData[tag], '__iter__'):
                            for data in sectionData[tag]:
                                dataString += "%s  %s\n" % (tag, data)
                        else:
                            dataString += "%s  %s\n" % (tag, sectionData[tag])
                        dataString += "XX\n"
            pwmData = None
            if pwm:
                pwmData = pwm
            else:
                pwmData = self.getMotifPWMClean(motif)
            if not pwmData: continue
            dataString += "P0  \tA\tC\tG\tT\n"
            for pwmIdx in range(len(pwmData)):
                pwmLine = "\t".join(map(str, pwmData[pwmIdx]))
                dataString += "%02i  \t%s\n" % (pwmIdx,pwmLine)
                
        return dataString                

    getSectionDataAsTransfac = getSectionData


    def exportDataAsPossum(self,motifname=None,species=None,pwm=None,AL="ACGT"):
        """ Export data to PoSSum format.
        For reference, see
        http://bibiserv.techfak.uni-bielefeld.de/possumsearch/files/PoSSuM-1_3.pdf

        =Sample BEGIN=
            BEGIN GROUP
            BEGIN INT
            ID Some matrix identifier
            AC Some accession
            DE A description describing the PSSM
            DE Multiple description lines are possible
            AP DNA
            LE 3
            # A C G T was specified by "AP DNA"
            MA 5 -1 -6 2
            MA -4 4 -1 -5
            MA 0 -3 3 -4
            END
    
            BEGIN FLOAT
            ID Some other matrix identifier
            DE Another description
            AL AUCG
            LE 2
            # A U C G
            MA 0.0 -3.5 3.2 -4.8
            MA -4.2 -1.0 4.0 -5.8
            END
    
            END
        =Sample END=

        Returns dataString as string type.

        >>> import motifparsers as M
        >>> d = M.parse_transfac("/home/prigor/projects/motifbrowser/rawdata/transfac/transfac_11.3.matrix.dat")
        >>> d.exportDataAsPossum("M00720")
        >>> print d.exportDataAsPossum("M00720")
        >>> outfile = file("/tmp/tst.file","w")
        >>> print >> outfile,d.exportDataAsPossum()
        """
        
        motifList = []
        if motifname==None:
            motifList = self.getMotifList(species)
        else:
            motifList.append(motifname)
        data = []
        
        #TODO: Need to determine FLOAT vs. INT values for the PWMs
        for motif in motifList:
            dataString = ""
            dataString += "ID %s\n" % motif
            dataString += "AC %s\n" % self.getMotifGeneName(motif)
            dataString += "AL %s\n" % AL
            pwmData = None

            if pwm:
                pwmData = pwm
            else:
                pwmData = self.getMotifPWMClean(motif)

            if not pwmData: continue

            dataString += "LE %s\n" % len(pwmData)
            for pwmIdx in range(len(pwmData)):
                pwmLine = " ".join(pwmData[pwmIdx])
                dataString += "MA %s\n" % pwmLine

            motifBlockPrepend = ""
            if "." in str(pwmData[0]):
                motifBlockPrepend = "BEGIN FLOAT"
            else:
                motifBlockPrepend = "BEGIN INT"

            dataString = "%s\n%s\nEND\n" % (motifBlockPrepend,dataString)
            data.append(dataString)

        #final enclosure
        dataString = "BEGIN GROUP\n%s\nEND\n" % "\n".join(data)
        return dataString                
                
    def getMotifGeneName(self,motifname):
        """ Need to overload this.
        
        """
        
        raise NotImplementedError
    
    def getMotifList(self,speciesName):
        """ Need to overload this.
        
        """
        
        raise NotImplementedError
    
    def getMotifPWMClean(self):
        """ Need to overload this.
        
        """
        
        raise NotImplementedError
    
    def getMotifConsensus(self,motifKey):
        """ Provides the consensus sequence using Xiaohui's code.
        Currently not in use because the stored matrices are counts, not frequencies or log_odds.
        """
        
        pwm = self.getMotifPWMClean(motifKey)
        sequence = calculateMotifConsensus(pwm)
        return sequence

############
### Parsers 
############


class TransfacParser(BaseMotifParser):
    """ Transfac database parser.
    """
    def __init__(self,filename,DEBUG=False):
        BaseMotifParser.__init__(self,filename,DEBUG)

    def generateData(self):
        """ Builds the motif dictionary.  Assumes no duplicate entries.
        """
        # State flags
        inData = False
        inMotif = False
        inSection = None
        inPWM = False

        # Special tags
        blankLineTag = "XX"
        sectionTag = "//"
        accTag = "AC"
        pwmTag = "P0"

        # Special holder
        motifName = None
        pwmKey = "pwm"  #used for the PWM key in order to normalize
                        #across different types of database sources

        # Let's do this...
        motifCount = 0
        sectionCount = 0
        for line in self.rawinputlines:
            tag = line[:2]
            data = str(line[3:]).strip() # remove leading and trailing spaces

            # Skip blanks, tagged or data
            if tag == blankLineTag: continue

            # Begin data area, always check for inData flag, otherwise skip line
            if tag == sectionTag and not inData:
                inData = True                
            if not inData: continue

            # Begin motif annotation
            if not inMotif and tag == sectionTag:
                sectionCount += 1
                inMotif = True
                if self.DEBUG:
                    print "In motif section: %s" % sectionCount
                continue

            # End motif annotation 
            if inMotif and tag==sectionTag:
                if self.DEBUG: 
                    print "Exiting motif: %s" % motifName
                #inMotif = False
                #motifName = None
                inPWM = False
                continue

            # Annotation section
            if inMotif and tag!=sectionTag:
                if tag==accTag and not inPWM:
                    motifName = data
                    self.motifs[motifName]={}
                    motifCount += 1
                    if self.DEBUG:
                        print "Entering motif: %s" % motifName
                        print "Total motifs found: %s" % motifCount
                    continue
                    
                if tag==pwmTag and not inPWM:
                    if self.DEBUG: print "\t\tEntering PWM of motif"
                    inPWM = True
                    self.motifs[motifName].update({pwmKey:[]})
                    continue

                # Cycle through PWM rows
                if inPWM:                    
                    try:
                        tag = int(tag)
                    except:
                        inPWM = False
                        if self.DEBUG: print "\t\tExiting PWM of motif"
                    if inPWM: 
                        if self.DEBUG: print "\t\tPWM:\t%s"%data
                        self.motifs[motifName][pwmKey].append(data)
                    continue
                
                # Handle other tags aside from Accession and PWM
                if not inPWM:
                    if data=="":
                        if self.DEBUG:
                            print "\t\tSkipping tag with blank data: %s" % tag
                        continue # skip blanks
                    if tag not in self.motifs[motifName].keys():
                        self.motifs[motifName].update({tag:[data]})
                    else:
                        self.motifs[motifName][tag].append(data)
                    if self.DEBUG: print "\t\t%s:\t%s"%(tag,data)

############
### Controllers
############

class TransfacController(BaseMotifController):
    """Controls transfac parser."""

    def __init__(self,filename,DEBUG=False):
        self.DEBUG = DEBUG
        self.parser = TransfacParser(filename,DEBUG)

    def parse(self):
        self.parser.loadFile()
        self.parser.generateData()

    def getData(self):
        return self.parser.getData()

class TransfacData(BaseMotifData):
    """Provides accessory functions to parsed data.
       Extend this class in order to add more 'tag' accessory methods.

       Currently supported tags: AC, NA, DE, BF, R?, P0 (pwm)
       
        In [14]: tfdata.getMotifList("yeast")
        Out[14]:
        ['M00170',
         'M00015',
         'M00031',
         'M00337',
         'M00213',
         'M00276',
         'M00274',
         'M00199',
         'M00197',
         'M01005',
         'M00288',
         'M00061',
         'M00168',
         'M00752',
         'M00049',
         'M00754',
         'M00125',
         'M01051',
         'M00303',
         'M00307',
         'M00713',
         'M00305',
         'M00207']
        
        In [15]: motif =  'M00276'
        
        In [16]: tfdata.getMotifReferences(motif)
        Out[16]: 'Cell differentiation by interaction of two HMG-box proteins: \
                 Mat1-Mc activates M cell-specific genes in S. pombe by recruiting \
                 the ubiquitous transcription factor Ste11 to weak binding sites\
                 PUBMED: 9233811.Kjaerulff S., Dooijes D., Clevers H., Nielsen O.EMBO J. 16:4021-4033 (1997).[1]; RE0006527.'
        
        In [17]: tfdata.getMotifDescription(motif)
        Out[17]: 'M-box interacting with Mat1-Mc'
        
        In [18]: tfdata.getMotifGeneName(motif)
        Out[18]: 'Mat1-Mc'
                
        In [19]: tfdata.getMotifPWM(motif)
        Out[19]:
        [['0', '4', '0', '6', 'Y'],
         ['0', '7', '0', '3', 'C'],
         ['2', '4', '1', '3', 'N'],
         ['9', '0', '0', '1', 'A'],
         ['1', '0', '0', '9', 'T'],
         ['0', '0', '0', '10', 'T'],
         ['0', '0', '10', '0', 'G'],
         ['0', '0', '0', '10', 'T'],
         ['1', '3', '0', '6', 'T'],
         ['3', '0', '2', '5', 'W']]
        
        In [20]: tfdata.getMotifPWMConsensus(motif)
        Out[20]: 'YCNATTGTTW'
        
        In [21]: tfdata.getMotifPWMClean(motif)
        Out[21]:
        [['0', '4', '0', '6'],
         ['0', '7', '0', '3'],
         ['2', '4', '1', '3'],
         ['9', '0', '0', '1'],
         ['1', '0', '0', '9'],
         ['0', '0', '0', '10'],
         ['0', '0', '10', '0'],
         ['0', '0', '0', '10'],
         ['1', '3', '0', '6'],
         ['3', '0', '2', '5']]
                  
    """
    
    def __init__(self,transfacdata,DEBUG=False):
        """Initializer.
        
        transfacdata contains the motif information returned by TransfacParser.getData()
        
        """
        
        BaseMotifData.__init__(self,transfacdata,DEBUG=DEBUG)
        
    def getMotifList(self,species=[]):
        """Get the motif list.

        Can also specify a species.

        """
        
        speciesTag = 'BF'
        motiflist = self.data.keys()
        if species:
            if type(species)==str:
                species=[species]
            species = [str(s).lower() for s in species]
            tempmotiflist = []
            for s in species:
                for motiflabel in motiflist:
                    try:
                        speciesSection = self.getSectionData(motiflabel,speciesTag)
                    except:
                        if self.DEBUG: print >>sys.stderr, "Excluding motif '%s' with no species information." % motiflabel
                        #tempmotiflist.append(motif)    
                        continue
                    for speciesEntry in speciesSection:
                        if s in str(speciesEntry).lower():
                            tempmotiflist.append(motiflabel)
            motiflist = self.uniquefy(tempmotiflist)
        return motiflist

    def getMotifSpeciesList(self,motiflabel,format=None):
        """Get the species list of both common and latin names. """
        def getSpecies(istring):
            """ Helper method to extract common and latin species names. """
            
            speclist = []
            
            res = istring.split("Species: ")
            
            if res and len(res)>1:
                res=res[1] # get the common,latin names
                if res[-1]==".":
                    res=res[:-1]
                speclist=res.split(",")
                speclist=[s.strip() for s in speclist]
            return speclist
            
        speciesTag = "BF"
        speciesSection = None
        try:
            speciesSection = self.getSectionData(motiflabel,speciesTag)
        except:
            if self.DEBUG: print >>sys.stderr, "Excluding motif '%s' with no species information." % motiflabel
        
        speclist = []
        if speciesSection:
            for speciesEntry in speciesSection:
                if speciesEntry:
                    speclist.append(getSpecies(speciesEntry))
            
        if format==str and speclist:
            speclist = ",".join(speclist)
        return speclist

    def getMotifGeneName(self,motifname):
        """Obtain associated gene name specified by the NA tag."""
        geneTag = "NA"
        genes = ""
        geneList = self.getSectionData(motifname,geneTag)
        if geneList: genes = ", ".join(geneList)       
        return genes            

    def getMotifDescription(self,motifname):
        """Obtain associated description specified by the DE tag."""
        descTag = "DE"
        descs = ""
        descEntries = self.getSectionData(motifname,descTag)
        if descEntries: descs = ", ".join(descEntries)
        return descs            

    def getMotifReferences(self,motifname):
        """Obtain references as a string using the transfac accession id, motifname."""
        self.is_valid_motif(motifname)
        refTags = ['RN','RA','RT','RL','RX']
        refs = ""
        for tag in self.getSectionList(motifname):
            if tag in refTags:
                refs += " ".join(self.getSectionData(motifname,tag))
        return refs                    

    def getMotifPWM(self,motifname):
        """Returns the pwm with default consensus vaules."""        
        pwmTag = "pwm"
        data = self.getSectionData(motifname,pwmTag)
        data = [line.split() for line in data]
        return data

    def getMotifPWMClean(self,motifname):
        """Returns the pwm matrix without default consensus values."""
        pwmKey = "pwm"  #used for the PWM key in order to normalize
                        #across different types of database sources

        data = self.getSectionData(motifname,pwmKey)
        data = [line.split()[:4] for line in data]
        return data

    def getMotifPWMLength(self,motifname):
        """Returns the pwm and length."""
        return len(self.getMotifPWMClean(motifname))

    def getMotifPWMConsensus(self,motifname):
        """Returns the default consensus sequences from transfac.
           Should overload this if a different consensus algorithm is desired.
        """
        pwm = self.getMotifPWM(motifname)
        seq = ""
        if str(pwm[0][-1]).isalpha():
            for index,entry in enumerate(pwm):
                if str(entry[-1]).isalpha():
                    seq+=entry[-1]
                else:
                    seq+="."
        else:
            seq = self.getMotifConsensus(motifname)
        return seq

# helper methods
def parse_transfac(filename,DEBUG=False):
    """Parses a transfac db file and returns a TransfacData instance.
    
    """
    
    tfcontrol = TransfacController(filename,DEBUG=DEBUG)
    tfcontrol.parse()
    return TransfacData(tfcontrol.getData(), DEBUG=DEBUG)

##################
# Parser registry
##################
dbparser = {}
dbparser["transfac"] = parse_transfac
