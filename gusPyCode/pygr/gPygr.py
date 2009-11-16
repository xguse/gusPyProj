"""
File contains my custom pygr classes and funcs.
"""

class solRead(object):
    '''INCOMPLETE'''
    recNum = 0
    def __init__(self,sortedFileLine,readLen=40):
        data = sortedFileLine.strip().split('\t')
        self.ID     = solRead.recNum+1
        solRead.recNum+=1
        
        self.name   = ':'.join(data[:7])
        self.seq    = data[8]
        self.qual   = data[9]
        self.chrm   = data[11]
        self.start  = data[12]
        self.end    = data[12]+readLen-1
        if data[13] == 'F':
            self.strand = 1
        elif data[13] == 'R':
            self.strand = -1
        else:
            print 'data[13] does not seem to contain R or F.  Your format may be wrong'
            exit()
