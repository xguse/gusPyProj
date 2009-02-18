import os
import sys
import pp
from time import time
from TAMO import MotifTools 
from TAMO.seq import Fasta 
from TAMO import MotifMetrics
from TAMO.MD.AlignAce import MetaAce 
from TAMO.MD.MDscan import MDscan 
from TAMO.MD.Meme import Meme 
import MarkovBackground
import seqStats

supportDirectory = 'path'
    
class AceWrap:
    """This is a wrapper to allow MDAP to initialize, run, and recieve results from TAMO's AlignAce
    interface.  It gets its optional aruments from those given to the MDAP main class by the user.
    The options are in the form an optparse options object (assume pass by ref).  This class should not be run 
    directly.  If you want to do this you should just call TAMO's interface directly.."""
    
    # Below is the __init__ def for TAMO's MetaAce class.
    # It is intended for help in maintaining this code.  However, you may need to look deeper
    # into TAMO.MD.AlignAce; I am afraid.
    # def __init__(self,fastafile = '', width=10, iterations = 5, gcback=0.38): 
    
    #-----------   
    def __init__(self, optionsObj, posArgs):
        """Sets the initial state of the wrapper class using user input from main script."""
        # Communicate user options to wrapper
        self.mdapOptions = optionsObj
        self.mdapArgs    = posArgs
        
        # object attributes
        self.output    = None
        self.dataStats = None
        
        #  Preliminary initialization of AlignAce specific options
        self.fastafile  = self.mdapArgs[1]
        self.width      = self.mdapOptions['kmerSize']
        self.iterations = self.mdapOptions['ace_iter']
        self.gcback     = None
    
    #-----------    
    def go(self):
        """Execution function: coordinates options used and background GC calculation, then runs
        TAMO.MD.AlignAce.MetaAce and catches the output in self.output for access from MDAP."""
        
        # Calc GC background of genomic sequences representing the
        # entire data set if requested.
        if self.mdapOptions['background'] == 1:
            self.dataStats = seqStats(self.mdapArgs[0])
            self.gcback = self.dataStats['percentGC']
            
        
        # call TAMO to do its thing
        self.output = MetaAce(self.fastafile, self.width, self.iterations, self.gcback)
        pass
    


#====================
class MemeWrap:
    """This is a wrapper to allow MDAP to initialize, run, and recieve results from TAMO's Meme
    interface.  It gets its optional aruments from those given to the MDAP main class by the user.
    The options are in the form an optparse options object (assume pass by ref).  This class should not be run 
    directly.  If you want to do this you should just call TAMO's interface directly."""
    
    # Below is the __init__ def for TAMO's Meme class:
    # It is intended for help in maintaining thos code.  However, you may need to look deeper
    # into TAMO.MD.Meme; I am afraid.
    # def __init__(self,file='', width='', extra_args='',genome='YEAST',bfile=None):
    
    #-----------   
    def __init__(self, optionsObj, posArgs):
        """Sets the initial state of the wrapper class using user input from main script."""
        
        # Communicate user options to wrapper
        self.mdapOptions = optionsObj
        self.mdapArgs    = posArgs
        
        # object attributes
        self.output    = None
        ##self.dataStats = None
        
        #  Preliminary initialization of Meme specific options
        self.file       = self.mdapArgs[1]  
        self.width      = self.mdapOptions['kmerSize'] # !--> will want to add support for ranged kmers.
        self.extra_args = None
        self.bfile      = None
        
        #-----------    
    def go(self):
        """Execution function: coordinates options used and any necessary background calculations,
        then runs TAMO.MD.Meme.Meme and catches the output in self.output for access from MDAP."""
        
        # if no bfile supplied call TAMO's markov background script and place results
        # in 'supportDirectory' with sensible name.
        if self.bfile == None:
            # Try to determine if there might already be a suitible bfile based on the name of the 
            # 'fastaOfAllSeqs' file supplied.  Otherwise create one and set bfileto point to it.
            if '%s.freq' % (self.mdapArgs[0]) in os.listdir(supportDirectory):
                sys.stderr.write("""ATTENTION: bfile found in %s that'fastaOfAllSeqs' file name.  This will be used. Check valididty of this file.""" \
                                 % (supportDirectory))
            else:
                MarkovBackground(self.mdapArgs[0],supportDirectory)
                sys.stderr.write("""ATTENTION: bfile <%s.freq> created and placed in %s.""" % \
                                 (self.mdapArgs[0].split('/')[-1], supportDirectory))
                self.bfile = '%s/%s.freq' % (supportDirectory, self.mdapArgs[0].split('/')[-1])
            
        
        # Call TAMO to do its thing:
        self.output = Meme(file=self.file, width=self.width, extra_args=self.extra_args, bfile=self.bfile)


#====================
class TamoWrap:
    """This is a wrapper to allow MDAP to initialize, run, and recieve results from TAMO's MotifMetrics
    exauhstive MD algorithm using its enrichment metric.  It gets its optional aruments from those given
    to the MDAP main class by the user.  The options are in the form an optparse options object
    (assume pass by ref).  This class should not be run directly.  If you want to do this, you should
    just call TAMO's interface directly."""
    
    # Below are the requirements for my script usingTAMO.MotifMetrics:
    # 1. fastaOfAllSeqs
    # 2. listOfLinkedSeqs
    # 3. kmerSize _OR_ kmerRange
    
    #-----------   
    def __init__(self, optionsObj, posArgs):
        """Sets the initial state of the wrapper class using user input from main script."""
        
        # Communicate user options to wrapper
        self.mdapOptions = optionsObj
        self.mdapArgs    = posArgs
        
        # object attributes
        self.output = None
        self.toFile = None
        
        #  Preliminary initialization of options
        self.allSeqs         = MotifMetrics.ProbeSet(self.mdapArgs[0])  
        self.linkedSeqs_ids  = map(lambda l: l.strip(), open(self.mdapArgs[1], 'rU')) # read geneNames into list killing \n's
        self.linkedSeqs_seqs = self.allSeqs.seqs_from_ids(self.linkedSeqs_ids) 
        self.kmerSize        = self.mdapOptions['kmerSize'] 
        if self.mdapOptions['kmerRange']:
            self.kmerRange =  map(lambda x: int(x), self.mdapOptions['kmerRange'].split(',')) # split the str(n,n) into a list of ints
        else:
            self.kmerRange = None
        
        #-----------    
    def go(self):
        """Execution function: coordinates options used then uses TAMO.MotifMetrics to
        find kmers with good enrichment in listOfLinkedSeqs. Catches the output in 
        self.output for access from MDAP."""
        
        # Initialize job_server for later
        job_server = pp.Server()
        
        # set metric thresholds here
        pVal_thresh     = 0.01
        church_thresh   = 0.01
        binomial_thresh = 0.01

        # # # # # # # # # # # # #
        # ::THIN THE HEARD PHASE::
        # Are we using a range or a single size? Then make a list of all kmers in range
        # that are present in at least 10% of linkedSeqs (top_nmers_seqs()) to reduce
        # needless kmer testing in the metrics phase.
        
        theShortList = []
        
        if self.kmerRange:
            for k in range(self.kmerRange[0],self.kmerRange[1]):
                kmers = MotifMetrics.top_nmers_seqs(k, self.linkedSeqs_seqs)
                print '%s %smers found.' % (len(kmers), k)
                theShortList.extend(kmers)
        else:
            theShortList = MotifMetrics.top_nmers_seqs(self.kmerSize, self.linkedSeqs_seqs)
            print '%s %smers found.' % (len(theShortList), self.kmerSize)
            
        # Convert theShortList into list of motif objs not just strings
        # REASON: church routine asks the motif for its width.
        for i in range(0,len(theShortList)):
            theShortList[i] = MotifTools.Motif_from_text(theShortList[i])
            
        # # # # # # # # # # # #
        # ::METRICS PHASE::
        # Using theShortList, calculate the:
        #       --------METRICS----------   --METHOD CALL--
        #     - HyperGeometric Enrichment      (p_value)
        #     - Group Specificity Score        (church)
        #     - Over-representation            (binomial)
        #
        # Retain those kmers that recieve the cut-off score or better in at least one
        # of the above metrics.
        
        # list with indexes as follows [kmer, p_value, church, binomial]
        keepers = []  
        
        t1 = time()
        count = 1
        shortList_Len = len(theShortList)
        for kmer in theShortList:
            p_value  = self.allSeqs.p_value(kmer, self.linkedSeqs_ids)
            church   = 'NA' #self.allSeqs.church(kmer, self.linkedSeqs_ids)
            binomial = 'NA' #self.allSeqs.binomial(kmer, self.linkedSeqs_ids)
            
            if p_value <= pVal_thresh or church <= church_thresh or binomial <= binomial_thresh:
                keepers.append([kmer, p_value, church, binomial])
                print '%s\t%s\t--\t%s of %s' % (kmer, p_value, count, shortList_Len)
            count+=1
        t2 = time()
        self.output = keepers
        print 'Calculating the metrics took %.3f min.' % ((t2-t1)/60) 
        
        # Create a formated string to be printed to a file in MDAP class.
        toFile = ['#kmer\tp_value\tchurch\tbinomial\n']
        for i in keepers:
            toFile.append('%s\t%s\t%s\t%s\n' % (i[0],i[1],i[2],i[3]))
            
        self.toFile = toFile
            