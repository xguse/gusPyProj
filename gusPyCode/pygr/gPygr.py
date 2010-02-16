"""
File contains my custom pygr classes and funcs.
"""
import pygr
from BCBio.GFF import GFFParser


class SliceData(object):
    '''ToDo:
    ___ enforce vaule types
    '''
    def __init__(self, seq_id, start, stop, orientation, xtrasDict={},**kwargs):
        self.seq_id      = seq_id
        self.start       = start
        self.stop        = stop
        self.orientation = orientation

        assert len(set(xtrasDict.keys()).intersection(set(self.__dict__.keys()))) == 0, \
               'ERROR: At least one key in xtrasDict clashes with existing obj attributes.'
        assert len(set(kwargs.keys()).intersection(set(self.__dict__.keys()))) == 0, \
               'ERROR: At least one key in kwargs clashes with existing obj attributes.'
        self.__dict__.update(xtrasDict)
        self.__dict__.update(kwargs)
        


class simpleGFF2PygrSQLite(object):
    def __init__(self,liteserver):
        """Uses BCBio.GFF.GFFparser.parse_simple() to initiate and populate SQLite tables.
        No 'nested' Parent:Child relationships are used, but P:C info can be accessed in parse_simple
        results.  For NOW, only essential data for pygrSliceObj is recorded.
        Returns: I THINK nothing --> external liteserver should be populated after loop?
        """
        self._parser   = GFFParser()
        self.srvrInfo  = liteserver
        self.tables    = {}  # <-- should add func to populate this on init if any already exist in serverInfo
        

    def _newRow(self, sliceObj):
        '''Creates new row in correct table.  Creates new table if need be.
        '''
        
        # Do we have to make a new table?
        if sliceObj.featureType in self.tables:
            # Get next k integer
            k = len(self.tables[sliceObj.featureType]) # works bc k starts w/ 0
            
            self.tables[sliceObj.featureType].new(k              = k,
                                                  seq_id         = sliceObj.seq_id,
                                                  start          = sliceObj.start,
                                                  stop           = sliceObj.stop,
                                                  orientation    = sliceObj.orientation,
                                                  featureType    = sliceObj.featureType,
                                                  gffKind        = sliceObj.gffKind,
                                                  coordReference = sliceObj.coordReference,
                                                  Parent         = sliceObj.Parent,
                                                  description    = sliceObj.description)
        else:
            self._newTable(sliceObj)
            
    def _newTable(self,sliceObj):
        '''Convert parsed GFF data to correct sql CREATE command, and add
        tableObj and keyName to self.tables.
        '''
        sqlCmd = 'CREATE TABLE %s (k INTEGER PRIMARY KEY, seq_id TEXT, start INT, stop INT, orientation INT, featureType TEXT, gffKind TEXT, coordReference TEXT, Parent TEXT, description TEXT,);'
        self.tables[sliceObj.featureType] = sqlgraph.SQLTable(sliceObj.featureType,
                                                              serverInfo=self.srvrInfo,
                                                              writeable=True,
                                                              createTable=sqlCmd)
        
        self.tables[sliceObj.featureType].new(k              = 0,
                                              seq_id         = sliceObj.seq_id,
                                              start          = sliceObj.start,
                                              stop           = sliceObj.stop,
                                              orientation    = sliceObj.orientation,
                                              featureType    = sliceObj.featureType,
                                              gffKind        = sliceObj.gffKind,
                                              coordReference = sliceObj.coordReference,
                                              Parent         = sliceObj.Parent,
                                              description    = sliceObj.description)
    def update(self,gffPath):
        gff      = open(gffFullPath)
        gffPrsr  = self._parser.parse_simple(gff)
        
        for rec in gffPrsr:
            # If rec == 'directive'  Move on.
            if rec.keys()[0] == 'directive':
                continue
            # else: create/add to correct table
            else:
                self.newRow(sGFFrow2SliceObj(gffRec)) # calls newTable if needed
        
        self.srvrInfo.commit()
        
    def getTableNames(self):
        srvrInfo.commit()
        cur = self.srvrInfo.cursor()
        x   = cur.execute('SELECT * FROM sqlite_master;')
        return [x[1] for x in x]    
    
    



##class solRead(object):
    ##'''INCOMPLETE'''
    ##recNum = 0
    ##def __init__(self,sortedFileLine,readLen=40):
        ##data = sortedFileLine.strip().split('\t')
        ##self.ID     = solRead.recNum+1
        ##solRead.recNum+=1
        
        ##self.name   = ':'.join(data[:7])
        ##self.seq    = data[8]
        ##self.qual   = data[9]
        ##self.chrm   = data[11]
        ##self.start  = data[12]
        ##self.end    = data[12]+readLen-1
        ##if data[13] == 'F':
            ##self.strand = 1
        ##elif data[13] == 'R':
            ##self.strand = -1
        ##else:
            ##print 'data[13] does not seem to contain R or F.  Your format may be wrong'
            ##exit()
            
# defs ##############
def sGFFrow2SliceObj(gffRec):
    """Takes a rec from BCBio.GFF.GFFParser.parse_simple() and returns a SliceData object.
    """
    recStart  = min(recData['location'])  # *
    recStop   = max(recData['location'])-1  # *
    recOri    = recData['strand']  # *
    recID     = recData['id']  # *
    
    recData   = gffRec[gffRec.keys()[0]][0]
    xtras={}
    xtras['gffKind']        = gffRec.keys()[0]
    xtras['featureType']    = recData['type']
    xtras['coordReference'] = recData['rec_id']
    
    if 'quals' in recData['quals']:
        if 'Parent' in recData['quals']:
            xtras['Parent'] = ','.join(recData['quals']['Parent'])
        else:
            xtras['Parent'] = None
        if 'description' in recData['quals']:
            xtras['description'] = ','.join(recData['quals']['description'])
        else:
            xtras['description'] = None
    else:
        xtras['Parent']      = None
        xtras['description'] = None
    
    
    
    return SliceData(recID,
                     recStart,
                     recStop,
                     recOri,
                     xtras)
